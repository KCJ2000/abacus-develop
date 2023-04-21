#include"../dnrm2.h"
#include"gtest/gtest.h"
#include"gmock/gmock.h"
#include<iostream>

/***************************************************
 *  unit test of double function dnrm2.
****************************************************/

/**
 * @brief Construct a new TEST object
 * -Tested fuction:double dnrm2(const int n, const double *x, const int incx) ;
 *      -int n is the maximum of the index
 *      -int incx is the interval of two index
 *      -the input is an double array type instead of std::vector type
 * The result is L2 norm
 */

TEST(dnrm2Test, dnrm2test){
    std::string break_string_out;
    testing::internal::CaptureStderr();
    double bg1[7]={1,2,3,4,5.3,5.4,9.24};
    int a1 = -1,b1 = -1;
    int a2 = 10,b2 = 1;
    double out;
    out = dnrm2(a1,bg1,b1);
    break_string_out = testing::internal::GetCapturedStderr();
    EXPECT_THAT(break_string_out,testing::HasSubstr(
        "\n error in dnrm2, n < 0 or incx <= 0, "));
    EXPECT_DOUBLE_EQ(out,0.0);

    out = dnrm2(a1,bg1,b2);
    break_string_out = testing::internal::GetCapturedStderr();
    EXPECT_THAT(break_string_out,testing::HasSubstr(
        "\n error in dnrm2, n < 0 or incx <= 0, "));
    EXPECT_DOUBLE_EQ(out,0.0);

    out = dnrm2(a2,bg1,b1);
    break_string_out = testing::internal::GetCapturedStderr();
    EXPECT_THAT(break_string_out,testing::HasSubstr(
        "\n error in dnrm2, n < 0 or incx <= 0, "));
    EXPECT_DOUBLE_EQ(out,0.0);

    out = dnrm2(a2,bg1,b2);
    break_string_out = testing::internal::GetCapturedStderr();
    EXPECT_THAT(break_string_out,testing::HasSubstr(
        "\n error in dnrm2, the index is out of the array"));
    EXPECT_DOUBLE_EQ(out,0.0);

    a2 = 0;
    b2 = 2;

    out = dnrm2(a2,bg1,b2);
    EXPECT_DOUBLE_EQ(out,0.0);

    a2 = 8;
    b2 = 2;
    out = dnrm2(a2,bg1,b2);
    double L2 = sqrt(bg1[0]*bg1[0]+bg1[2]*bg1[2]+
                        bg1[4]*bg1[4]+bg1[6]*bg1[6]);
    EXPECT_DOUBLE_EQ(out , L2);

}
