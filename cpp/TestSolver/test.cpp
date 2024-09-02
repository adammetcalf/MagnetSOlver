#include "pch.h"																				//unit test libraries
#include <gtest/gtest.h>																		//unit test libraries
#include <iostream>

#include "../MagnetSolver/Joint.h"																//include header file for joint class
#include "../MagnetSolver/Tentacle.h"															//include header file for tentacle class


// Define the test fixture for TestJoint
class TestJoint : public ::testing::Test {
protected:
    
    Joint* joint;                                                                               // Null Pointer to a Joint object

    void SetUp() override {                                                                     // SETUP ---- Runs before Each test
        
        joint = new Joint(0.0f, 90.0f, 0.0f, 1.0f);                                            // Initialize the joint object with some test values
    }
 
    void TearDown() override {                                                                  // TEARDOWN --- Runs after each test
        
        delete joint;                                                                           // Delete the joint object
    }
};

// Define the test fixture for TestTentacle
class TestTentacle : public ::testing::Test {
protected:

    Tentacle* tentacle;                                                                         // null pointer to tentacle object

    void SetUp() override {

        Eigen::Matrix<float, 7, 2> Angles;
        Angles << 90.0, 180.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0,
            0.0, 0.0;

        tentacle = new Tentacle(0.01f, Angles);
    };

    void TearDown() override {

        delete tentacle;
    }
};

// Test to verify that the joint object is constructed properly
TEST_F(TestJoint, TestConstruction) {
    EXPECT_TRUE(joint != nullptr);                                                              // Check that the joint exists
}

// Test to check the correctness of the DH frame after construction
TEST_F(TestJoint, TestFrameAfterConstruction) {
    Eigen::Matrix4f frame = joint->getFrame();                                                  // Get the joint frame

    Eigen::Matrix4f expAns;                                                                     // Define Matehmatically expected Matrix
    expAns << 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 1.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0;

    //std::cout << expAns << std::endl;                                                         // display the expected answer
    //std::cout << frame << std::endl;                                                          // display the actual answer
    EXPECT_EQ(frame, expAns);
}

TEST_F(TestJoint, TestUpdateAngles) {
    Eigen::Matrix4f oldframe = joint->getFrame();                                               // Get the joint frame

    Eigen::Matrix4f expAns;                                                                     // Define Matehmatically expected Matrix
    expAns << 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, -1.0, 0.0,
        0.0, 1.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0;

    EXPECT_EQ(oldframe, expAns);                                                                // Check correct init

    joint->UpdateAngles({ 90.0f,90.0f });                                                       // call update angles, with [90,90]

    Eigen::Matrix4f newframe = joint->getFrame();                                               // Get the new joint frame

    EXPECT_NE(oldframe, newframe);
    EXPECT_NE(newframe, expAns);

    Eigen::Matrix4f NewexpAns;                                                                  // Define Matehmatically expected Matrix for updated angles
    NewexpAns << 0.0, 0.0, 1.0, 0.0,
        1.0, 0.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0,
        0.0, 0.0, 0.0, 1.0;

    EXPECT_EQ(newframe, NewexpAns);
}

TEST_F(TestTentacle, TestConstruction) {

    EXPECT_TRUE(tentacle != nullptr);                                                           // Check that the tentacle exists

    Eigen::Matrix<float, Eigen::Dynamic, 2> JointAngles = tentacle->GetAngles();                // get the current angles

    //std::cout << JointAngles << std::endl;                                                    // display the joint angles

    Eigen::Matrix<float, 7, 2> NewAngles;
    NewAngles << 90.0, 180.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0,
        0.0, 0.0;

    EXPECT_EQ(JointAngles, NewAngles);                                                          // Check NewANgles == JointAngles
}

TEST_F(TestTentacle, TestUpdateAngles) {

    Eigen::Matrix<float, Eigen::Dynamic, 2> JointAngles = tentacle->GetAngles();                // Get the current angles

    Eigen::Matrix<float, 7, 2> NewAngles;
    NewAngles << 90.0, 180.0,
        0.0, 90.0,
        0.0, 90.0,
        0.0, 90.0,
        0.0, 90.0,
        0.0, 90.0,
        0.0, 90.0;

    tentacle->UpdateAngles(NewAngles);                                                          // Update the Angles

    Eigen::Matrix<float, Eigen::Dynamic, 2> NewJointAngles = tentacle->GetAngles();             // Get the updated angles

    EXPECT_EQ(NewJointAngles, NewAngles);
}