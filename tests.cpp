#define BOOST_TEST_MODULE CPITests
#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_CASE(test_addition) {
    int result = 2 + 3;
    BOOST_CHECK_EQUAL(result, 5);
}