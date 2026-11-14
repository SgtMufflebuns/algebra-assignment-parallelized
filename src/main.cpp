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

//Function name could possibly change if we introduce PLU Decomposition. 
Matrix cyclicDecomposition(Matrix matrix, int rank, int size){
    //Implement cyclic decomposition using rank

    //Here is where we will use the information of rank and size to help compute determinant

    //const auto matrixBlockSize = dimension / (size - 1); This does not work if processes do not neatly match up with  

    if(rank > 0){
        //Do something here
    }

    return matrix;
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

    // Initialize MPI environment
    //Initialize MPI without argc and argv since we are not accepting any arguments.
    MPI_Init(nullptr, nullptr);  

    //Determine the rank and the size of the specific process.
    int rank, size; // every rank needs to know who it is before deciding what work to do
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // root is used to determine if the current process is the main process responsible for output and coordination
    const bool root = rank == 0;

        
    try {
        const auto total_start = Clock::now();
        Field::initialize();
        if (root) {
            std::cout << std::fixed << std::setprecision(6)
                      << "Matrix dimension " << dimension << '\n'
                      << "Seed " << seed << '\n'
                      << "Field order " << Field::cardinality() << '\n'
                      << "Modulus " << Field::modulus() << '\n';
        }

        if (root) {
            std::cout << "Generating A and B\n" << std::flush;
        }
        auto start = Clock::now();
        Matrix a;
        Matrix b;
        populate(a, b);
        if (root) {
            std::cout << "Generation seconds " << seconds(start) << '\n';
        }

        if (root) {
            std::cout << "Computing C = A * B\n" << std::flush;
        }
        start = Clock::now();
        const Matrix c = a.multiply(b);
        if (root) {
            std::cout << "Product seconds " << seconds(start) << '\n';
        }



        cyclicDecomposition(a, rank, size);
        cyclicDecomposition(b, rank, size);

        if (root) {
            std::cout << "Computing det(A)\n" << std::flush;
        }
        start = Clock::now();
        const auto det_a = determinant(a);
        if (root) {
            print_determinant("det(A)", det_a, seconds(start));
        }

        if (root) {
            std::cout << "Computing det(B)\n" << std::flush;
        }
        start = Clock::now();
        const auto det_b = determinant(b);
        if (root) {
            print_determinant("det(B)", det_b, seconds(start));
        }

        if (root) {
            std::cout << "Computing det(C)\n" << std::flush;
        }
        start = Clock::now();
        const auto det_c = determinant(c);
        if (root) {
            print_determinant("det(C)", det_c, seconds(start));
        }

        

        const auto expected = det_a.multiply(det_b);
        const bool verified = det_c == expected;
        if (root) {
            std::cout << "det(A) * det(B) = " << expected << '\n'
                      << "Identity verified " << std::boolalpha << verified << '\n'
                      << "Total seconds " << seconds(total_start) << '\n';
        }

        MPI_Finalize(); // should be the last MPI call a rank makes
        return verified ? 0 : 1;
    } catch (const std::exception& error) {
        std::cerr << "rank " << rank << ": " << error.what() << '\n';
        MPI_Abort(MPI_COMM_WORLD, 1);
        return 1;
    }
}
