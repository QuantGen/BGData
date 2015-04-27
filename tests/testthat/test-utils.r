library(BGData)

context("utils")

test_that("normalizeType", {
    
    expect_equal(typeof(normalizeType('double')), 'double')
    expect_equal(typeof(normalizeType(double())), 'double')
    expect_equal(typeof(normalizeType('integer')), 'integer')
    expect_equal(typeof(normalizeType(integer())), 'integer')
    expect_equal(typeof(normalizeType('character')), 'character')
    expect_equal(typeof(normalizeType(character())), 'character')
    expect_equal(typeof(normalizeType('complex')), 'complex')
    expect_equal(typeof(normalizeType(complex())), 'complex')
    expect_equal(typeof(normalizeType('test')), 'character')
    expect_equal(typeof(normalizeType(1)), 'double')
    expect_equal(typeof(normalizeType(1L)), 'integer')
    
})
