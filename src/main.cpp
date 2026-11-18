#include <mpi.h>
#include <chrono>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

#include <algebra/finite_fields/finite_fields.hpp>
#include <algebra/matrices/matrices.hpp>

namespace {

using Field = algebra::FiniteField<11, 3>;
using Element = Field::Element;
constexpr std::size_t dimension = ALGEBRA_ASSIGNMENT_SIZE;
constexpr std::uint32_t seed = 123456u;
using Matrix = algebra::Matrix<Element, dimension, dimension>;
using Clock = std::chrono::steady_clock;

static_assert(algebra::FieldElement<Element>);
static_assert(algebra::RingElement<Matrix>);

double seconds(Clock::time_point start) {
    return std::chrono::duration<double>(Clock::now() - start).count();
}

std::uint32_t draw_coefficient(std::mt19937& generator, std::uint32_t prime) {
    if (prime < 2) {
        throw std::invalid_argument("The prime must be at least two");
    }
    constexpr std::uint64_t range = std::uint64_t{1} << 32;
    const std::uint64_t limit = range - range % prime;
    for (;;) {
        const std::uint64_t value = generator();
        if (value < limit) {
            return static_cast<std::uint32_t>(value % prime);
        }
    }
}

Element next_element(std::mt19937& generator) {
    Element::polynomial_type polynomial;
    for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        const auto coefficient = draw_coefficient(generator, Field::characteristic);
        polynomial.set_coefficient(degree, Element::coefficient_type::from_integer(coefficient));
    }
    return Element::from_polynomial(polynomial);
}

void populate(Matrix& a, Matrix& b) {
    std::mt19937 generator(seed);
    for (std::size_t i = 0; i < dimension; ++i) {
        for (std::size_t j = 0; j < dimension; ++j) {
            a(i, j) = next_element(generator);
            b(i, j) = next_element(generator);
        }
    }
}

int calculateOwner(int globalColumn, int numNodes)
{
    // calculate the owner of a specific column,
    // in this instance, for the 1D column cyclic,
    // is easily done with modulo operations

    return globalColumn % numNodes;
}

template<class Element>
void broadcastFiniteFieldElement(Element& element, int owner, MPI_Comm commWorld) 
{
    // Definitely the trickiest part here; problem is that
    // element can't cleanly be broadcast as a single MPI datatype.
    // So we need to break it into its coefficients and broadcast each one.

    // need this in order to extract the # of degrees
    using Field = typename Element::parent_type;

    // these are the actual values of the element
    using Coefficient = typename Element::coefficient_type;
    using Polynomial = typename Element::polynomial_type;

    int rank;
    MPI_Comm_rank(commWorld, &rank);

    constexpr std::size_t numDegrees = Field::extension_degree;
    
    //create a buffer for the coefficients
    //to be stored in for the network...
    std::array<std::uint32_t, numDegrees> coefficientBuffer{};

    if (rank == owner) {
        for (std::size_t i = 0; i < numDegrees; ++i) {
            coefficientBuffer[i] = element.coefficient(i).representative();
        }
    }

    //broadcast this buffer to all nodes...
    MPI_Bcast(
        coefficientBuffer.data(),
        static_cast<int>(numDegrees),
        MPI_UINT32_T,
        owner,
        commWorld
    );

    //...and reconstruct the element using the buffer.
    if (rank != owner) 
    {
        // construct a temp polynomial to replace the new element. 
        Polynomial tempPolynomial;
        // since polynomials are dynamic, can increase size with no issue.
        for (std::size_t i = 0; i < numDegrees; ++i) 
        {
            // i'th degree with value of the i'th coefficient in the buffer
            tempPolynomial.set_coefficient(i, Coefficient::from_index(coefficientBuffer[i]));
        }
        element = Element::from_polynomial(tempPolynomial);
    }
}

void broadcastFiniteFieldMultipliers(std::vector<Element>& multipliers, int owner, MPI_Comm commWorld)
{
    // Implementation for broadcasting multipliers en masse
    
    // need this in order to extract the # of degrees
    using Field = typename Element::parent_type;

    // these are the actual values of the element
    using Coefficient = typename Element::coefficient_type;
    using Polynomial = typename Element::polynomial_type;

    int rank;
    MPI_Comm_rank(commWorld, &rank);

    constexpr std::size_t numDegrees = Field::extension_degree;
    
    // create a larger buffer to store the coefficients
    // for every single multiplier element to batch
    // send in one go
    std::vector<std::uint32_t> allCoefficientsBuffer(multipliers.size() * numDegrees, 0);
    
    if(rank == owner)
    {
        for(std::size_t i = 0; i < multipliers.size(); ++i)
        {
            for(std::size_t j = 0; j < numDegrees; ++j)
            {
                allCoefficientsBuffer[i * numDegrees + j] = multipliers[i].coefficient(j).representative();
            }
        }
    }

    //broadcast this buffer to all nodes...
    MPI_Bcast(
        allCoefficientsBuffer.data(),
        static_cast<int>(multipliers.size() * numDegrees),
        MPI_UINT32_T,
        owner,
        commWorld
    );

    //...and reconstruct the multipliers using the buffer.
    if(rank != owner)
    {
        for(std::size_t i = 0; i < multipliers.size(); ++i)
        {
            // construct a temp polynomial to replace the new element. 
            Polynomial tempPolynomial;
            // since polynomials are dynamic, can increase size with no issue.
            for(std::size_t j = 0; j < numDegrees; ++j)
            {
                // now need to skip over the coefficients of the previous multipliers, so we get the right corresponding one
                tempPolynomial.set_coefficient(j, Coefficient::from_index(allCoefficientsBuffer[i * numDegrees + j]));
            }
            multipliers[i] = Element::from_polynomial(tempPolynomial);
        }
    }

}

Element determinant(Matrix matrix, MPI_Comm commWorld) {

    // we will be parallelizing this using 1D-column cyclic, 
    // which is more efficient compared to the blocked implementation.

    int rank, size;
    MPI_Comm_rank(commWorld, &rank);
    MPI_Comm_size(commWorld, &size);

    // with 1D-column cyclic, the owner of a column alternates
    // between each node to help with load balancing.

    /* ex: 2 nodes, 4x4 matrix 

        0 1 0 1
        0 1 0 1
        0 1 0 1
        0 1 0 1 
    */

    Element result = Element::one();
    for (std::size_t activeColumn = 0; activeColumn < dimension; ++activeColumn) {
        std::size_t pivot_row = -1;
        std::size_t owner = calculateOwner(activeColumn, size);

        // owner of the column handle the pivot row
        if(rank == owner)
        {
            pivot_row = activeColumn;
            while (pivot_row < dimension && matrix(pivot_row, activeColumn).is_zero()) 
            {
                ++pivot_row;
            }
        }

        // broadcast the pivot row to all nodes
        MPI_Bcast(&pivot_row, 1, MPI_UNSIGNED_LONG, owner, commWorld);

        if (pivot_row == dimension) 
        {
            return Element::zero();
        }

        if (pivot_row != activeColumn) 
        {
            matrix.swap_rows(pivot_row, activeColumn);
            result = result.negate();
        }

        // since this pivot can be stale on other columns,
        // broadcast it from the owner (its accurate)
        Element pivot = Element::zero();
        Element inverse = Element::one();
        int remainingRows = dimension - activeColumn - 1;

        // need to also create a vector for the multipliers
        // for the owner to track and distribute
        std::vector<Element> multipliers(remainingRows, Element::zero());
        
        // the owner of the pivot column grabs the 
        // true pivot value and calculates the multipliers
        if (rank == owner) 
        {
            pivot = matrix(activeColumn, activeColumn);
            if (!pivot.is_one()) 
            {
                inverse = pivot.reciprocal();
            }

            for (std::size_t targetRow = activeColumn + 1; targetRow < dimension; ++targetRow) 
            {
                if (!matrix(targetRow, activeColumn).is_zero()) 
                {
                    multipliers[targetRow - activeColumn - 1] =
                        matrix(targetRow, activeColumn).multiply(inverse);
                }
            }
        }

        // now that the pivot & multipliers have been set, 
        // need to broadcast it to all nodes by reconstructing
        // each finite field element and rebuilding piecemeal
        broadcastFiniteFieldElement(pivot, owner, commWorld);
        broadcastFiniteFieldMultipliers(multipliers, owner, commWorld);

        result = result.multiply(pivot);

        // finally, paralellization takes place; each node handles
        // multiplying ONLY the columns assigned to them.
        // this means each node only has to deal with 1/n of this work
        // hopefully, this offsets the costs of MPI
        for (std::size_t targetCol = rank; targetCol < dimension; targetCol += size) 
        {
            if (targetCol <= activeColumn) 
            {
                continue;
            }

            // remember that we're making an upper triangular matrix, 
            // so we only need to worry about the elements above the diagonal
            for (std::size_t targetRow = activeColumn + 1; targetRow < dimension; ++targetRow) 
            {
                const Element multiplier = multipliers[targetRow - activeColumn - 1];

                if (multipliers[targetRow - activeColumn - 1].is_zero()) 
                {
                    continue;
                }

                const Element correction = multiplier.multiply(matrix(activeColumn, targetCol));
                matrix(targetRow, targetCol) = matrix(targetRow, targetCol).add(correction.negate());
            }
        }
        
        
    }
    return result;
}

void print_determinant(const char* name, const Element& value, double elapsed, int rank) {
    if(rank == 0)
    {
        std::cout << name << " = " << value << '\n';
        std::cout << "Coefficients [";
        for (std::size_t degree = 0; degree < Field::extension_degree; ++degree) {
        if (degree != 0) {
            std::cout << ", ";
        }
        std::cout << value.coefficient(degree);
    }
    std::cout << "]\nSeconds " << elapsed << '\n' << std::flush;
    }
}

}

double CalculateSlowestElapsedTime(double localElapsedTime, MPI_Comm commWorld) 
{
    double globalElapsedTime = 0.0f;
    MPI_Allreduce(&localElapsedTime, &globalElapsedTime, 1, MPI_DOUBLE, MPI_MAX, commWorld);

    return globalElapsedTime;
}

int main(int argc, char* argv[]) {
    bool mpi_initialized = false;
    try {
        if (MPI_Init(&argc, &argv) != MPI_SUCCESS) 
        {
            throw std::runtime_error("MPI_Init failed!");
        }
        mpi_initialized = true;
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        // everything's passed, can get to work

        double total_start = MPI_Wtime();

        char host[MPI_MAX_PROCESSOR_NAME] = {};
        int length;
        MPI_Get_processor_name(host, &length);

        Field::initialize();

        if (rank == 0) 
        {
            std::cout << std::fixed << std::setprecision(6)
                      << "Matrix dimension " << dimension << '\n'
                      << "Seed " << seed << '\n'
                      << "Field order " << Field::cardinality() << '\n'
                      << "Modulus " << Field::modulus() << '\n';

            std::cout << "Generating A and B\n" << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();
        Matrix a;
        Matrix b;
        populate(a, b);
        double localElapsed = MPI_Wtime() - start;
        double globalElapsed = CalculateSlowestElapsedTime(localElapsed, MPI_COMM_WORLD);
        if (rank == 0) 
        {
            std::cout << "Generation seconds " << globalElapsed << '\n';
            std::cout << "Computing C = A * B\n" << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const Matrix c = a.multiply(b);
        localElapsed = MPI_Wtime() - start;
        globalElapsed = CalculateSlowestElapsedTime(localElapsed, MPI_COMM_WORLD);
        if (rank == 0) 
        {
            std::cout << "Product seconds " << globalElapsed << '\n';
            std::cout << "Computing det(A)\n" << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const auto det_a = determinant(a, MPI_COMM_WORLD);
        localElapsed = MPI_Wtime() - start;
        globalElapsed = CalculateSlowestElapsedTime(localElapsed, MPI_COMM_WORLD);
        print_determinant("det(A)", det_a, globalElapsed, rank);

        if (rank == 0) 
        {
            std::cout << "Computing det(B)\n" << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const auto det_b = determinant(b, MPI_COMM_WORLD);
        localElapsed = MPI_Wtime() - start;
        globalElapsed = CalculateSlowestElapsedTime(localElapsed, MPI_COMM_WORLD);
        print_determinant("det(B)", det_b, globalElapsed, rank);

        if (rank == 0) 
        {
            std::cout << "Computing det(C)\n" << std::flush;
        }

        MPI_Barrier(MPI_COMM_WORLD);
        start = MPI_Wtime();
        const auto det_c = determinant(c, MPI_COMM_WORLD);
        localElapsed = MPI_Wtime() - start;
        globalElapsed = CalculateSlowestElapsedTime(localElapsed, MPI_COMM_WORLD);
        print_determinant("det(C)", det_c, globalElapsed, rank);

        const auto expected = det_a.multiply(det_b);
        const bool verified = det_c == expected;
        if (rank == 0) {
            std::cout << "det(A) * det(B) = " << expected << '\n'
                      << "Identity verified " << std::boolalpha << verified << '\n'
                      << "Total seconds " << MPI_Wtime() - total_start << '\n';
        }
        if (MPI_Finalize() != MPI_SUCCESS) 
        {
            throw std::runtime_error("MPI_Finalize failed");
        }
        mpi_initialized = false;
        return verified ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        if (mpi_initialized) {
            MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        }

        return 1;
    }
    return 0;
}
