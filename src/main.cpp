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
#include <mpi.h>

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

void serialize(
    Element element,
    std::vector<uint32_t> &v
){
    //Take in an element and then return a uint32_t array as the output
    std::array<algebra::PrimeFieldElement<11>, 3> arr; 
    arr = element.coefficients();

    for(algebra::PrimeFieldElement<11> pElement: arr){
        std::uint32_t value = pElement.representative();
        v.push_back(value);
    }
}

void cyclicBlockDistribution(
    std::size_t blockSize, 
    int(&coords) [2], 
    int (&processSpace) [2],
    int& processRows,
    int& processCols
){
    
    int localRows = 0;
    int localCols = 0;

    //Get the number of rows and cols each chunk requires. 
    //Combines cyclical and block distribution.
    for(std::size_t i = 0; i < dimension; ++i){
        int blockRow = i / blockSize; 

        if(blockRow % processRows == coords[0]){
            ++localRows;
        }
    }

    for(std::size_t j = 0; j < dimension; ++j){
        int blockCol = j / blockSize;

        if(blockCol % processCols == coords[1]){
            ++localCols;
        }
    }

    processSpace[0] = localRows;
    processSpace[1] = localCols; 
}

//Function name could possibly change if we introduce PLU Decomposition. 
void cyclicDecompositionDistribution(Matrix matrix, int size){
    //Implement cyclic decomposition using rank

    //Here is where we will use the information of rank and size to help compute determinant
    MPI_Comm cart_comm; 

    //Arrange the processes P in a Pr x Pc grid.
    const int dimSize = 2;
    int dims[dimSize]{0, 0};
    int periods[dimSize]{0, 0};
    int coords[dimSize]{};
    MPI_Dims_create(size, dimSize, dims);

    int processRows = dims[0]; // Pr
    int processCols = dims[1]; // Pc
    int reorder = 0; 

    MPI_Cart_create(MPI_COMM_WORLD, dimSize, dims, periods, reorder, &cart_comm);

    int cartRank;
    MPI_Comm_rank(cart_comm, &cartRank);
    MPI_Cart_coords(cart_comm, cartRank, dimSize, coords);

    const auto blockSize = 4; //Number of elements in each block MB or NB
    
    //We need to allocate memory for each local chunk 
    std::size_t numBlocks = (dimension + blockSize - 1) / blockSize; 

    int processSpace[2]{};

    //Get the number of Rows and Columns for the speciifc processes' 
    cyclicBlockDistribution(blockSize, processSpace, coords, processRows, processCols);

    int localRows = processSpace[0];
    int localCols = processSpace[1];

    std::vector<Element> processMatrix(localRows * localCols); 

    //Processing block by block
    for(std::size_t block = 0; block < numBlocks * numBlocks; ++block){
        
        //Row major order calc for block position. 
        std::size_t blockRow = block / numBlocks;
        std::size_t blockCol = block % numBlocks; 
        
        //We can do this cause of MPI_Cart calcs earlier
        int processRowCoordinate = blockRow % processRows;
        int processColCoordinate = blockCol % processCols; 

        //Accounts for uneven blocks where the matrix may not evenly divide. 
        std::size_t rowEnd = std::min((blockRow + 1) * blockSize, dimension);
        std::size_t colEnd = std::min((blockCol + 1) * blockSize, dimension);

        int position = 0;
        std::size_t messageSize = sizeof(std::uint32_t) * 3 * blockSize * blockSize;
        char buff[messageSize];

        //Loops through each element in the block
        for(std::size_t i = blockRow * blockSize; i < rowEnd; ++i){
            for(std::size_t j = blockCol * blockSize; j < colEnd;  ++j){
                //Pass in an element at a time first.
                std::vector<std::uint32_t> coeffs; 
                
                //Serialize the element before passing
                serialize(matrix(i, j), coeffs);

                MPI_Pack(&coeffs, 3, MPI_INT, buff, messageSize, &position, cart_comm);
            }
        }
        
        int processCartRank;
        int processCoords[2]{};
        processCoords[0] = processRowCoordinate;
        processCoords[1] = processColCoordinate; 

        MPI_Cart_rank(cart_comm, processCoords, &processCartRank);

        //If the rank is yourself then dont send the data. 
        if(processCartRank == 0){
            MPI_Send(buff, position, MPI_PACKED, processCartRank, 0, cart_comm);
        }

        /* if(processRowCoordinate == coords[0] && processColCoordinate == coords[1]){
            
            //What row and column in terms of blocks is this block located. 
            int localBlockRow = blockRow / processRows;
            int localBlockCol = blockCol / processCols;

            //Iterate over each element in the block
            for(std::size_t globalRow = blockRow * blockSize; globalRow < rowEnd; ++globalRow){
                for(std::size_t globalCol = blockCol * blockSize; globalCol < colEnd; ++globalCol){
                    std::size_t localRow = localBlockRow * blockSize + globalRow % blockSize;
                    std::size_t localCol = localBlockCol * blockSize + globalCol % blockSize; 

                    processMatrix[localRow * localCols + localCol] = matrix(globalRow, globalCol);
                }
            }

        } */
        
    }
}

Element determinant(Matrix matrix) {
    Element result = Element::one();
    for (std::size_t k = 0; k < dimension; ++k) {
        std::size_t pivot_row = k;
        while (pivot_row < dimension && matrix(pivot_row, k).is_zero()) {
            ++pivot_row;
        }
        if (pivot_row == dimension) {
            return Element::zero();
        }
        if (pivot_row != k) {
            matrix.swap_rows(pivot_row, k);
            result = result.negate();
        }
        const Element pivot = matrix(k, k);
        result = result.multiply(pivot);
        Element inverse = Element::one();
        if (!pivot.is_one()) {
            inverse = pivot.reciprocal();
        }
        for (std::size_t i = k + 1; i < dimension; ++i) {
            if (matrix(i, k).is_zero()) {
                continue;
            }
            const Element multiplier = matrix(i, k).multiply(inverse);
            matrix(i, k) = Element::zero();
            for (std::size_t j = k + 1; j < dimension; ++j) {
                const Element correction = multiplier.multiply(matrix(k, j));
                matrix(i, j) = matrix(i, j).add(correction.negate());
            }
        }
    }
    return result;
}

void print_determinant(const char* name, const Element& value, double elapsed) {
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

int main() {
    try {
        const auto total_start = Clock::now();
        Field::initialize();
        std::cout << std::fixed << std::setprecision(6)
                  << "Matrix dimension " << dimension << '\n'
                  << "Seed " << seed << '\n'
                  << "Field order " << Field::cardinality() << '\n'
                  << "Modulus " << Field::modulus() << '\n';

        std::cout << "Generating A and B\n" << std::flush;
        auto start = Clock::now();
        Matrix a;
        Matrix b;
        populate(a, b);
        std::cout << "Generation seconds " << seconds(start) << '\n';

        std::cout << "Computing C = A * B\n" << std::flush;
        start = Clock::now();
        const Matrix c = a.multiply(b);
        std::cout << "Product seconds " << seconds(start) << '\n';

        //Initialize MPI without argc and argv since we are not accepting any arguments.
        MPI_Init(nullptr, nullptr);  

        //Determine the rank and the size of the specific process.
        int rank, size;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

        if(rank == 0){
            cyclicDecompositionDistribution(a, size);
        }

        std::cout << "Computing det(A)\n" << std::flush;
        start = Clock::now();
        const auto det_a = determinant(a);
        print_determinant("det(A)", det_a, seconds(start));

        std::cout << "Computing det(B)\n" << std::flush;
        start = Clock::now();
        const auto det_b = determinant(b);
        print_determinant("det(B)", det_b, seconds(start));

        std::cout << "Computing det(C)\n" << std::flush;
        start = Clock::now();
        const auto det_c = determinant(c);
        print_determinant("det(C)", det_c, seconds(start));

        MPI_Finalize();

        const auto expected = det_a.multiply(det_b);
        const bool verified = det_c == expected;
        std::cout << "det(A) * det(B) = " << expected << '\n'
                  << "Identity verified " << std::boolalpha << verified << '\n'
                  << "Total seconds " << seconds(total_start) << '\n';
        return verified ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
