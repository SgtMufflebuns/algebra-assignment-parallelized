include(FetchContent)

if(NOT TARGET algebra::finite_fields)
    FetchContent_Declare(
        algebra_finite_fields
        GIT_REPOSITORY git@github.com:CalvinWu2020/algebra-finite-fields-group-1.git
        GIT_TAG 4cf1e522f3774a2820b16b0a38a4f07cfe12686e
    )
    FetchContent_MakeAvailable(algebra_finite_fields)
endif()

if(NOT TARGET algebra::matrices)
    FetchContent_Declare(
        algebra_matrices
        GIT_REPOSITORY git@github.com:nathanjrussell/algebra-matrices.git
        GIT_TAG 84abddac26d2805072efcb1c44bbe7faac45fa6b
    )
    FetchContent_MakeAvailable(algebra_matrices)
endif()
