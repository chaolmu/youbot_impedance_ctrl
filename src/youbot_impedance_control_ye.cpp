//============================================================================
// Name        : youbot_impedance_control_ye.cpp
// Author      : tianyu
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include<iostream>
#include <eigen3/Eigen/Dense>
#include <math.h>
#include <fstream>
#include <ros/ros.h>
#include <brics_actuator/JointPositions.h>
#include <brics_actuator/JointTorques.h>
#include <brics_actuator/JointVelocities.h>
#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
#include <sensor_msgs/JointState.h>
#include <impedance_control_sub.h>

#include <dynamic_reconfigure/server.h>
#include <youbot_impedance_control_ye/ImpedanzeParaConfig.h>


using namespace Eigen;

//const Vector2d a(0.135, 0.218);
//const double pi = 3.141592653589793;
//const double joint_offsets[5] = {2.9496, 1.1344, -2.5481, 1.7889, 2.9234};
//const Vector2d X_soll(0.15,0.2);

Vector2d NullTorque(0, 0);
const double zeitschritt = 0.003;//second

ros::Publisher armPositionsPublisher;
ros::Publisher armTorquePublisher;
ros::Publisher armVelocitiesPublisher;
ros::Subscriber joint_statesSubscriber;

ros::Publisher armTorqueBeobachter;
int main(int argc, char** argv)
{
	ros::init(argc, argv, "youbot_impedanz_controll_ye");
	ros::NodeHandle nh;

	armPositionsPublisher = nh.advertise<brics_actuator::JointPositions>("arm_1/arm_controller/position_command", 1);
	armTorquePublisher = nh.advertise<brics_actuator::JointTorques>("arm_1/arm_controller/torques_command", 1);
	armTorqueBeobachter =nh.advertise<brics_actuator::JointTorques>("Torque_Beobachter", 1);

	dynamic_reconfigure::Server<youbot_impedance_control_ye::ImpedanzeParaConfig> server;
	dynamic_reconfigure::Server<youbot_impedance_control_ye::ImpedanzeParaConfig>::CallbackType f;
	f = boost::bind(&Impecallback, _1, _2);
	server.setCallback(f);

	ros::Rate rate(1/zeitschritt);

	//subscribe jointstate actual_p actual_v
	Vector2d READY_sub(0, 0);
	Subscri sub(READY_sub);
	joint_statesSubscriber = nh.subscribe("joint_states", 1, &Subscri::callback,&sub);

	Vector2d theta_y(-2, 0.05);
	Vector2d theta_y_dot(1, -1);
	Vector2d theta_c(0, 0);
	Vector2d theta_c_dot(0, 0);
	Vector2d theta_I(0, 0);
	Vector2d theta_I_dot(0, 0);

	Vector2d X_actual(0,0);
	Vector2d tau(0,0);
	Vector2d y(0,0);
	Vector2d pub_torques(0,0);
	Matrix<double,2,2> Ja;
	Matrix<double,2,2> Ja_dot;
	Matrix<double,2,2> Ja_inverse;
	Ja.setZero();

	int wait=0;
	std::ofstream XY_Schwingung;
	XY_Schwingung.open("XY_Schwingung.txt", std::ios::out);

	while (nh.ok() && wait < (1/zeitschritt))  //1s
	{
		wait = wait + 1;
		//	ros::spinOnce();
		rate.sleep();
	}

	while(nh.ok()&&Stop==false)
	{
		theta_y=sub.actual_p;
		theta_y_dot=sub.actual_v;
//		std::cout<<"	theta_y="<<theta_y.transpose();
		youbot2torque(theta_y,theta_y_dot,theta_c,theta_c_dot);

		youbot2Impedanz(theta_y,theta_y_dot,theta_I,theta_I_dot);

				DirekteKinematik(theta_I,X_actual);
//				std::cout<<"	X_actual="<<X_actual.transpose();

				getJa(theta_I,Ja);
				getJa_dot(theta_I,theta_I_dot,Ja_dot);
				getJa_Inverse(Ja,Ja_inverse);
				ImpedanzRegler(theta_I,theta_I_dot,X_actual,Ja,Ja_dot,Ja_inverse,y);
//				std::cout<<"	Y="<<y.transpose();
				Inverse_Dynamic(theta_c,theta_c_dot,y,tau);

				// if tau in intervall else max or min
				tau=tau.array().min(tau_max.array());
				tau=tau.array().max(tau_min.array()); 
				pub_torques=tau;//aufpassen! vorzeichen in MsgGenerator geaendert.

				armTorquePublisher.publish(generateJointTorqueMsg(pub_torques));
				armTorqueBeobachter.publish(generateJointTorqueMsg(X_actual));

//		std::cout<<'\n';
		ros::spinOnce();
		rate.sleep();

		writeToFile(XY_Schwingung,X_actual);

	}
	armTorquePublisher.publish(generateJointTorqueMsg(NullTorque));
	XY_Schwingung.close();
	std::cout<<"finish to end!!";
	return 0;
}
