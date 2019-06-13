context("utils")

test_that("normalizeType", {

    expect_equal(typeof(BGData:::normalizeType("double")), "double")
    expect_equal(typeof(BGData:::normalizeType(double())), "double")
    expect_equal(typeof(BGData:::normalizeType("integer")), "integer")
    expect_equal(typeof(BGData:::normalizeType(integer())), "integer")
    expect_equal(typeof(BGData:::normalizeType("character")), "character")
    expect_equal(typeof(BGData:::normalizeType(character())), "character")
    expect_equal(typeof(BGData:::normalizeType("complex")), "complex")
    expect_equal(typeof(BGData:::normalizeType(complex())), "complex")
    expect_warning(BGData:::normalizeType("test"))
    expect_equal(suppressWarnings(typeof(BGData:::normalizeType("test"))), "character")
    expect_equal(typeof(BGData:::normalizeType(1)), "double")
    expect_equal(typeof(BGData:::normalizeType(1L)), "integer")

})
