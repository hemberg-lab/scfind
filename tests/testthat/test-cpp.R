test_that("testing eliasFanoCoding()", {
    res <- eliasFanoCoding(
        list(c(5, 8, 8, 15, 32), c(2, 3, 5, 7, 11, 13, 24)), 
        c(2,2)
    )
    expect_is(res, "list")
    expect_is(res$H, "list")
    expect_is(res$L, "list")
    expect_equal(length(res), 2)
    expect_equal(as.numeric(res$H[[1]]), c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1))
    expect_equal(as.numeric(res$L[[1]]), c(0, 1, 0, 0, 0, 0, 1, 1, 0, 0))
    expect_equal(as.numeric(res$H[[2]]), c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1))
    expect_equal(as.numeric(res$L[[2]]), c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0))
    expect_error(res <- eliasFanoCoding(list(a, b)))
    expect_error(res <- eliasFanoCoding("test", l))
})

test_that("testing eliasFanoDecoding()", {
    H1 <- c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1)
    L1 <- c(0, 1, 0, 0, 0, 0, 1, 1, 0, 0)
    H2 <- c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1)
    L2 <- c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0)
    res <- eliasFanoDecoding(
        c(0, 1, 0, 1, 1, 0, 1, 0, 0, 0, 0, 0, 1), 
        c(0, 1, 0, 0, 0, 0, 1, 1, 0, 0), 
        2
    )
    expect_is(res, "numeric")
    expect_equal(length(res), 5)
    expect_equal(res, c(5, 8, 8, 15, 32))
    
    res <- eliasFanoDecoding(
        c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1), 
        c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0), 
        2
    )
    expect_is(res, "numeric")
    expect_equal(length(res), 7)
    expect_equal(res, c(2, 3, 5, 7, 11, 13, 24))
    
    expect_error(res <- eliasFanoDecoding(
        c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1))
    )
    expect_error(res <- eliasFanoDecoding(
        c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1),
        c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0)
    ))
    expect_error(res <- eliasFanoDecoding(
        yan,
        c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
        2
    ))
    expect_error(res <- eliasFanoDecoding(
        c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1),
        yan,
        2
    ))
    expect_error(res <- eliasFanoDecoding(
        c(1, 1, 0, 1, 1, 0, 1, 0, 1, 0, 0, 0, 1),
        c(1, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0),
        yan
    ))
})

test_that("testing int2bin()", {
    expect_is(int2bin(2), "logical")
    expect_equal(as.numeric(int2bin(2)), c(1, 0))
    expect_equal(as.numeric(int2bin(1)), 1)
    expect_error(int2bin("test"))
})
