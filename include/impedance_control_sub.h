#ifndef _IMPEDANCE_CONTROL_SUB_H_
#define _IMPEDANCE_CONTROL_SUB_H_

#include<iostream>
#include <eigen3/Eigen/Dense>
#include <values_internet.h>
#include <math.h>
#include <fstream>
#include <ros/ros.h>
#include <brics_actuator/JointPositions.h>
#include <brics_actuator/JointTorques.h>
#include <brics_actuator/JointVelocities.h>
#include <boost/units/systems/si.hpp>
#include <boost/units/io.hpp>
#include <sensor_msgs/JointState.h>

#include <dynamic_reconfigure/server.h>
#include <youbot_impedance_control_ye/ImpedanzeParaConfig.h>


using namespace Eigen;

typedef Matrix<double,5,1>  Vector5d;
const Vector2d a(0.135, 0.218);
const double pi = 3.141592653589793;
const double joint_offsets[5] = {2.9496, 1.1344, -2.5481, 1.7889, 2.9234};

//const Vector2d X_soll(0.15,0.2);
Vector2d X_soll(0.12,0.25);
Vector2d K_p(120,120);
Vector2d K_d(10,10);
Vector2d M_d(0.5,0.5);

const Vector2d Xdot_soll(0,0);
const Vector2d Xddot_soll(0,0);
const Vector2d tau_max(6.0,2);
const Vector2d tau_min(-6.0,-2);
bool Stop=false;

void Impecallback(youbot_impedance_control_ye::ImpedanzeParaConfig &config, uint32_t level) 
{
			X_soll<< config.X,   config.Y;
			K_p<<config.Kp_X, config.Kp_Y;
			K_d<<config.Kd_X, config.Kd_Y;
			M_d<<config.Md_X, config.Md_Y;
			M_d=M_d.cwiseInverse();
//			std::cout<<"	M_d="<<M_d.transpose();
//			std::cout<<'\n';
			//M_d=M_d.cwiseInverse();
			Stop=config.Stop;
}

void getJa(Vector2d &theta, Matrix<double,2,2> &Ja)
{
	double s1 = sin(theta(0));
	double s12 = sin(theta(0)+ (theta(1)));
	double c1 = cos(theta(0));
	double c12 = cos(theta(0) + (theta(1)));

	Ja << -a(0) * s1 - a(1) * s12, -a(1) * s12,
			a(0) * c1 + a(1) * c12, a(1) * c12;
}
void getJa_dot(Vector2d &theta, Vector2d &theta_dot, Matrix<double,2,2> &Ja_dot)
{
		double s1 = sin(theta(0));
		double c1 = cos(theta(0));
		double s12 = sin(theta(0) + (theta(1)));
		double c12 = cos(theta(0) + (theta(1)));

		Ja_dot << -a(0) * c1 * theta_dot(0) - a(1) * c12 * (theta_dot(0) + (theta_dot(1))), -a(1) * c12 * (theta_dot(0) + (theta_dot(1))),
				   -a(0) * s1 * theta_dot(0) - a(1) * s12 * (theta_dot(0) + (theta_dot(1))), -a(1) * s12 * (theta_dot(0) + (theta_dot(1)));
}
void getJa_Inverse(Matrix<double,2,2> Ja, Matrix<double,2,2> &Ja_inverse)
{
	Ja_inverse=Ja.inverse();
}
void DirekteKinematik(Vector2d &theta, Vector2d &X_actual)
{
	X_actual << a(0) * cos(theta(0) ) + a(1) * cos(theta(0)  + (theta(1))),
		 a(0) * sin(theta(0) ) + a(1) * sin(theta(0)  + (theta(1)));
}
void ImpedanzRegler(Vector2d &theta_I,Vector2d &theta_I_dot, Vector2d &X_actual,
		Matrix<double,2,2> &Ja,Matrix<double,2,2> &Ja_dot,Matrix<double,2,2> &Ja_inverse,Vector2d &y)
{
	Vector2d X_dot=Ja*theta_I_dot;	
	Vector2d X_ddot=Ja_dot*theta_I_dot;

	y=(((X_soll.array()-X_actual.array())*K_p.array()+(Xdot_soll.array()-X_dot.array())*K_d.array())*M_d.array()+(Xddot_soll.array()-X_ddot.array()));

	y=Ja_inverse*y;//output from Controller
    y(0)=-y(0);//input for calculation
}
void getM(Vector2d & theta, Matrix<double,2,2> & M)
{
	  
  double sin00 = sin((-2.9234));
  double cos00 = cos((-2.9234));
  double cos01 = cos(theta(1));
  double cos02 = cos((-1.1344) - theta(0) + theta(1) - (-2.9234));
  double sin01 = sin((-1.1344) - theta(0) + theta(1) + (-2.9234));
  double cos03 = cos((-1.1344) - theta(0) + theta(1) + (-2.9234));
  double sin02 = sin((-1.1344) - theta(0) + theta(1) - (-2.9234));
  double cos04 = cos((-1.1344) - theta(0) + theta(1) - (2.0 * (-2.9234)));
  double cos05 = cos((-1.1344) - theta(0) + theta(1) + (2.0 * (-2.9234)));
  double cos06 = cos(theta(0) - theta(1));
  double sin03 = sin(theta(1));
  double cos07 = cos((-1.1344) - theta(0) + theta(1));
  double sin04 = sin(theta(0) - theta(1));
  double sin05 = sin((-1.1344) - theta(0) - (-2.9234));
  double cos08 = cos((-1.1344) - theta(0) + (-2.9234));
  double cos09 = cos((-1.1344) - theta(0) - (-2.9234));
  double sin06 = sin((-1.1344) - theta(0) + (-2.9234));
  double sin07 = sin((-1.1344) - theta(0) + theta(1) - (2.0 * (-2.9234)));
  double sin08 = sin((-1.1344) - theta(0) + theta(1) + (2.0 * (-2.9234)));
  double cos10 = cos((-1.1344) - theta(0));
  double cos11 = cos((-1.1344) + (-2.9234));
  double sin09 = sin((-1.1344) - (-2.9234));
  double sin10 = sin((-1.1344) + (-2.9234));
  double cos12 = cos((-1.1344) - (-2.9234));
  double sin11 = sin((-1.1344) - theta(0) + theta(1));
  double cos13 = cos((-1.1344));
  double sin12 = sin((-1.1344) - theta(0));
  double sin13 = sin((-1.1344));
  double sin14 = sin(theta(0));
  double cos14 = cos(theta(0));	  
  
	M(0,0) = I3y + I4y + (pow(c3x,2.0) * l3m) + (pow(c3z,2.0) * l3m) + (l4m * pow(c4z + (l4z * cos01),2.0)) + (l4m * pow(c4x + (l4z * sin03),2.0)) + (I5y * pow(cos00,2.0)) + (Igy * pow(cos00,2.0)) + (l5m * pow(l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03),2.0)) + (mg * pow(l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03),2.0)) + (I5x * pow(sin00,2.0)) + (Igx * pow(sin00,2.0)) + (l5m * pow(cos00,2.0) * pow(c5z + l5z + (l4z * cos01),2.0)) + (mg * pow(cos00,2.0) * pow(cgz + l5z + (l4z * cos01),2.0)) + (l5m * pow(sin00,2.0) * pow(c5z + l5z + (l4z * cos01),2.0)) + (mg * pow(sin00,2.0) * pow(cgz + l5z + (l4z * cos01),2.0));
	M(0,1) = -I4y - (I5y * pow(cos00,2.0)) - (Igy * pow(cos00,2.0)) - (I5x * pow(sin00,2.0)) - (Igx * pow(sin00,2.0)) - (c4z * l4m * (c4z + (l4z * cos01))) - (l5m * (l5x + (c5x * cos00) - (c5y * sin00)) * (l5x + (c5x * cos00) - (c5y * sin00) + (l4z * sin03))) - (mg * (l5x + (cgx * cos00) - (cgy * sin00)) * (l5x + (cgx * cos00) - (cgy * sin00) + (l4z * sin03))) - (c4x * l4m * (c4x + (l4z * sin03))) - (l5m * pow(cos00,2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01))) - (mg * pow(cos00,2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01))) - (l5m * pow(sin00,2.0) * (c5z + l5z) * (c5z + l5z + (l4z * cos01))) - (mg * pow(sin00,2.0) * (cgz + l5z) * (cgz + l5z + (l4z * cos01)));
	M(1,0) = M(0,1);
	M(1,1) = I4y + (l5m * pow(l5x + (c5x * cos00) - (c5y * sin00),2.0)) + (mg * pow(l5x + (cgx * cos00) - (cgy * sin00),2.0)) + (pow(c4x,2.0) * l4m) + (pow(c4z,2.0) * l4m) + (I5y * pow(cos00,2.0)) + (Igy * pow(cos00,2.0)) + (I5x * pow(sin00,2.0)) + (Igx * pow(sin00,2.0)) + (l5m * pow(cos00,2.0) * pow(c5z + l5z,2.0)) + (mg * pow(cos00,2.0) * pow(cgz + l5z,2.0)) + (l5m * pow(sin00,2.0) * pow(c5z + l5z,2.0)) + (mg * pow(sin00,2.0) * pow(cgz + l5z,2.0));


}
void getC(Vector2d & theta,Vector2d & theta_dot, Matrix<double,2,2> & C)
{
	  
  double sin00 = sin(2.0 * (-2.9234));
  double sin01 = sin((-1.1344) - theta(0) + theta(1) - (2.0 * (-2.9234)));
  double sin02 = sin((-1.1344) - theta(0) + theta(1) + (2.0 * (-2.9234)));
  double cos00 = cos((-1.1344) - theta(0) + theta(1) - (-2.9234));
  double sin03 = sin((-1.1344) - theta(0) + theta(1) + (-2.9234));
  double cos01 = cos((-1.1344) - theta(0) + theta(1) + (-2.9234));
  double sin04 = sin((-1.1344) - theta(0) + theta(1) - (-2.9234));
  double sin05 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)));
  double sin06 = sin((-1.1344) - theta(0) + theta(1));
  double cos02 = cos((-2.9234));
  double sin07 = sin((-2.9234));
  double sin08 = sin(theta(1));
  double sin09 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) + (2.0 * (-2.9234)));
  double sin10 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) - (2.0 * (-2.9234)));
  double cos03 = cos(theta(1));
  double sin11 = sin(theta(0) - theta(1));
  double sin12 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) - (-2.9234));
  double sin13 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) + (-2.9234));
  double cos04 = cos((-1.1344) - theta(0) + theta(1));
  double cos05 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) - (-2.9234));
  double cos06 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) + (-2.9234));
  double cos07 = cos(theta(1) + (-2.9234));
  double cos08 = cos((-1.1344) - theta(0) + theta(1) - (2.0 * (-2.9234)));
  double cos09 = cos(theta(1) - (-2.9234));
  double sin14 = sin(theta(1) + (-2.9234));
  double sin15 = sin(theta(1) - (-2.9234));
  double cos10 = cos((-1.1344) - theta(0) + theta(1) + (2.0 * (-2.9234)));
  double cos11 = cos(2.0 * (-2.9234));
  double cos12 = cos(theta(0) - theta(1));
  double sin16 = sin(theta(0) - theta(1) + (-2.9234));
  double cos13 = cos(theta(0) - theta(1) - (-2.9234));
  double cos14 = cos(theta(0) - theta(1) + (-2.9234));
  double sin17 = sin(theta(0) - theta(1) - (-2.9234));
  double sin18 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)));
  double sin19 = sin((2.0 * (-1.1344)) - theta(0) + theta(1));
  double cos15 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)));
  double sin20 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1));
  double cos16 = cos((-1.1344) - theta(0) + (-2.9234));
  double sin21 = sin((-1.1344) - theta(0) + (-2.9234));
  double sin22 = sin((-1.1344) - theta(0) - (-2.9234));
  double cos17 = cos((-1.1344) - theta(0) - (-2.9234));
  double cos18 = cos((-1.1344) - theta(0));
  double sin23 = sin(2.0 * (-1.1344));
  double sin24 = sin(theta(0));
  double cos19 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1));
  double cos20 = cos((2.0 * (-1.1344)) - theta(0) + theta(1));
  double cos21 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1) - (-2.9234));
  double sin25 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1) - (-2.9234));
  double sin26 = sin((2.0 * (-1.1344)) - theta(0) + theta(1) + (-2.9234));
  double cos22 = cos((2.0 * (-1.1344)) - theta(0) + theta(1) - (-2.9234));
  double cos23 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) - (2.0 * (-2.9234)));
  double sin27 = sin((2.0 * (-1.1344)) - theta(0) + theta(1) - (-2.9234));
  double sin28 = sin((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1) + (-2.9234));
  double cos24 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + theta(1) + (-2.9234));
  double cos25 = cos((2.0 * (-1.1344)) - theta(0) + theta(1) + (-2.9234));
  double cos26 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)) + (2.0 * theta(1)) + (2.0 * (-2.9234)));
  double sin29 = sin((2.0 * (-1.1344)) - theta(0));
  double sin30 = sin((-1.1344) - theta(0));
  double sin31 = sin((-1.1344) + (-2.9234));
  double cos27 = cos((-1.1344));
  double cos28 = cos((-1.1344) + (-2.9234));
  double sin32 = sin((-1.1344) - (-2.9234));
  double cos29 = cos((-1.1344) - (-2.9234));
  double sin33 = sin((-1.1344));
  double cos30 = cos(theta(0));
  double cos31 = cos((2.0 * (-1.1344)) - (2.0 * theta(0)));
  double cos32 = cos((2.0 * (-1.1344)) - theta(0));
  double cos33 = cos(2.0 * (-1.1344));
	  
	C(0,0) = -(0 * ((I5y * sin00) - (I5x * sin00) - (Igx * sin00) + (Igy * sin00) + (2.0 * l5m * ((c5y * cos02) + (c5x * sin07)) * (l5x + (c5x * cos02) - (c5y * sin07) + (l4z * sin08))) + (2.0 * mg * ((cgy * cos02) + (cgx * sin07)) * (l5x + (cgx * cos02) - (cgy * sin07) + (l4z * sin08)))) / 2.0) - (l4z * theta_dot(1) * ((c4z * l4m * sin08) - (l5m * l5x * cos03) - (l5x * mg * cos03) - (c4x * l4m * cos03) + (c5z * l5m * sin08) + (cgz * mg * sin08) + (l5m * l5z * sin08) + (l5z * mg * sin08) - (c5x * l5m * cos03 * cos02) - (cgx * mg * cos03 * cos02) + (c5y * l5m * cos03 * sin07) + (cgy * mg * cos03 * sin07)));
	C(0,1) = (I5y * 0 * sin00 / 2.0) - (I5x * 0 * sin00 / 2.0) - (Igx * 0 * sin00 / 2.0) + (Igy * 0 * sin00 / 2.0) - (c5x * c5y * l5m * 0) - (cgx * cgy * mg * 0) + (pow(c5x,2.0) * l5m * 0 * sin00 / 2.0) - (pow(c5y,2.0) * l5m * 0 * sin00 / 2.0) + (pow(cgx,2.0) * mg * 0 * sin00 / 2.0) - (pow(cgy,2.0) * mg * 0 * sin00 / 2.0) + (l4z * l5z * mg * 0 * sin08) - (l4z * l5z * mg * theta_dot(0) * sin08) + (l4z * l5z * mg * theta_dot(1) * sin08) + (2.0 * c5x * c5y * l5m * 0 * pow(cos02,2.0)) + (2.0 * cgx * cgy * mg * 0 * pow(cos02,2.0)) + (c5y * l5m * l5x * 0 * cos02) - (c4x * l4m * l4z * 0 * cos03) + (c4x * l4m * l4z * theta_dot(0) * cos03) - (c4x * l4m * l4z * theta_dot(1) * cos03) + (cgy * l5x * mg * 0 * cos02) - (l5m * l5x * l4z * 0 * cos03) + (l5m * l5x * l4z * theta_dot(0) * cos03) - (l5m * l5x * l4z * theta_dot(1) * cos03) - (l5x * l4z * mg * 0 * cos03) + (l5x * l4z * mg * theta_dot(0) * cos03) - (l5x * l4z * mg * theta_dot(1) * cos03) + (c5x * l5m * l5x * 0 * sin07) + (c4z * l4m * l4z * 0 * sin08) - (c4z * l4m * l4z * theta_dot(0) * sin08) + (c4z * l4m * l4z * theta_dot(1) * sin08) + (c5z * l5m * l4z * 0 * sin08) - (c5z * l5m * l4z * theta_dot(0) * sin08) + (c5z * l5m * l4z * theta_dot(1) * sin08) + (cgx * l5x * mg * 0 * sin07) + (cgz * l4z * mg * 0 * sin08) - (cgz * l4z * mg * theta_dot(0) * sin08) + (cgz * l4z * mg * theta_dot(1) * sin08) + (l5m * l4z * l5z * 0 * sin08) - (l5m * l4z * l5z * theta_dot(0) * sin08) + (l5m * l4z * l5z * theta_dot(1) * sin08) - (c5x * l5m * l4z * 0 * cos03 * cos02) + (c5x * l5m * l4z * theta_dot(0) * cos03 * cos02) - (c5x * l5m * l4z * theta_dot(1) * cos03 * cos02) - (cgx * l4z * mg * 0 * cos03 * cos02) + (cgx * l4z * mg * theta_dot(0) * cos03 * cos02) - (cgx * l4z * mg * theta_dot(1) * cos03 * cos02) + (c5y * l5m * l4z * 0 * cos03 * sin07) - (c5y * l5m * l4z * theta_dot(0) * cos03 * sin07) + (c5y * l5m * l4z * theta_dot(1) * cos03 * sin07) + (c5y * l5m * l4z * 0 * cos02 * sin08) + (cgy * l4z * mg * 0 * cos03 * sin07) - (cgy * l4z * mg * theta_dot(0) * cos03 * sin07) + (cgy * l4z * mg * theta_dot(1) * cos03 * sin07) + (cgy * l4z * mg * 0 * cos02 * sin08) + (c5x * l5m * l4z * 0 * sin08 * sin07) + (cgx * l4z * mg * 0 * sin08 * sin07);
	C(1,0) = (I5y * 0 * sin00 / 2.0) - (I5x * 0 * sin00 / 2.0) - (Igx * 0 * sin00 / 2.0) + (Igy * 0 * sin00 / 2.0) - (c5x * c5y * l5m * 0) - (cgx * cgy * mg * 0) + (pow(c5x,2.0) * l5m * 0 * sin00 / 2.0) - (pow(c5y,2.0) * l5m * 0 * sin00 / 2.0) + (pow(cgx,2.0) * mg * 0 * sin00 / 2.0) - (pow(cgy,2.0) * mg * 0 * sin00 / 2.0) - (l4z * l5z * mg * 0 * sin08) + (l4z * l5z * mg * theta_dot(0) * sin08) + (2.0 * c5x * c5y * l5m * 0 * pow(cos02,2.0)) + (2.0 * cgx * cgy * mg * 0 * pow(cos02,2.0)) + (c5y * l5m * l5x * 0 * cos02) + (c4x * l4m * l4z * 0 * cos03) - (c4x * l4m * l4z * theta_dot(0) * cos03) + (cgy * l5x * mg * 0 * cos02) + (l5m * l5x * l4z * 0 * cos03) - (l5m * l5x * l4z * theta_dot(0) * cos03) + (l5x * l4z * mg * 0 * cos03) - (l5x * l4z * mg * theta_dot(0) * cos03) + (c5x * l5m * l5x * 0 * sin07) - (c4z * l4m * l4z * 0 * sin08) + (c4z * l4m * l4z * theta_dot(0) * sin08) - (c5z * l5m * l4z * 0 * sin08) + (c5z * l5m * l4z * theta_dot(0) * sin08) + (cgx * l5x * mg * 0 * sin07) - (cgz * l4z * mg * 0 * sin08) + (cgz * l4z * mg * theta_dot(0) * sin08) - (l5m * l4z * l5z * 0 * sin08) + (l5m * l4z * l5z * theta_dot(0) * sin08) + (c5x * l5m * l4z * 0 * cos03 * cos02) - (c5x * l5m * l4z * theta_dot(0) * cos03 * cos02) + (cgx * l4z * mg * 0 * cos03 * cos02) - (cgx * l4z * mg * theta_dot(0) * cos03 * cos02) - (c5y * l5m * l4z * 0 * cos03 * sin07) + (c5y * l5m * l4z * theta_dot(0) * cos03 * sin07) - (cgy * l4z * mg * 0 * cos03 * sin07) + (cgy * l4z * mg * theta_dot(0) * cos03 * sin07);
	C(1,1) = -(0 * ((I5y * sin00) - (I5x * sin00) - (Igx * sin00) + (Igy * sin00) + (2.0 * l5m * ((c5y * cos02) + (c5x * sin07)) * (l5x + (c5x * cos02) - (c5y * sin07))) + (2.0 * mg * ((cgy * cos02) + (cgx * sin07)) * (l5x + (cgx * cos02) - (cgy * sin07)))) / 2.0);

}
void getN(Vector2d & theta, Vector2d &N)
{

	double sin00 = sin(-1.1344 - theta(0) + theta(1));
	double cos00 = cos(-1.1344 - theta(0) + theta(1));
	double sin01 = sin(-1.1344 - theta(0));
	double cos01 = cos(-1.1344 - theta(0) + theta(1) +2.9234);
	double cos02 = cos(-1.1344 - theta(0) + theta(1) -2.9234);
	double sin02 = sin(-1.1344 - theta(0) + theta(1) +2.9234);
	double sin03 = sin(-1.1344 - theta(0) + theta(1) -2.9234);
	double sin04 = sin(-1.1344);
//N(0)*0.67
	N(0) = (981.0 * c5y * l5m * sin03 / 200.0) - (981.0 * cgx * mg * cos02 / 200.0) - (981.0 * c5x * l5m * cos02 / 200.0) + (981.0 * cgy * mg * sin03 / 200.0) - (981.0 * c3x * l3m * cos(theta(0) - theta(1)) / 100.0) + (981.0 * c3z * l3m * sin01 / 100.0) + (981.0 * l4m * l4z * sin01 / 100.0) + (981.0 * l5m * l4z * sin01 / 100.0) + (981.0 * l4z * mg * sin01 / 100.0) - (981.0 * c5x * l5m * cos01 / 200.0) - (981.0 * cgx * mg * cos01 / 200.0) - (981.0 * c5y * l5m * sin02 / 200.0) - (981.0 * cgy * mg * sin02 / 200.0) - (981.0 * c4x * l4m * cos00 / 100.0) - (981.0 * l5m * l5x * cos00 / 100.0) - (981.0 * l5x * mg * cos00 / 100.0) + (981.0 * c4z * l4m * sin00 / 100.0) + (981.0 * c5z * l5m * sin00 / 100.0) + (981.0 * cgz * mg * sin00 / 100.0) + (981.0 * l5m * l5z * sin00 / 100.0) + (981.0 * l5z * mg * sin00 / 100.0);
	
	N(1) = (981.0 * c5x * l5m * cos02 / 200.0) + (981.0 * cgx * mg * cos02 / 200.0) - (981.0 * c5y * l5m * sin03 / 200.0) - (981.0 * cgy * mg * sin03 / 200.0) + (981.0 * c5x * l5m * cos01 / 200.0) + (981.0 * cgx * mg * cos01 / 200.0) + (981.0 * c5y * l5m * sin02 / 200.0) + (981.0 * cgy * mg * sin02 / 200.0) + (981.0 * c4x * l4m * cos00 / 100.0) + (981.0 * l5m * l5x * cos00 / 100.0) + (981.0 * l5x * mg * cos00 / 100.0) - (981.0 * c4z * l4m * sin00 / 100.0) - (981.0 * c5z * l5m * sin00 / 100.0) - (981.0 * cgz * mg * sin00 / 100.0) - (981.0 * l5m * l5z * sin00 / 100.0) - (981.0 * l5z * mg * sin00 / 100.0);

	//N(1)=N(1)*1.3;
}
void Inverse_Dynamic(Vector2d &theta,Vector2d &theta_dot,Vector2d &y,Vector2d &tau)
{
	Matrix<double,2,2> M;
	Matrix<double,2,2> C;
	Vector2d g_v(0,0);
	Vector2d F_v(0,0);
	getM(theta,M);
    getC(theta,theta_dot,C);
    getN(theta,g_v);

		//if(theta(0)<(-2.03))
		if(theta_dot(0)>0)
		{
		F_v(0)=0.4;
		}
		//else if(theta(0)>(-1.98))
		else if(theta_dot(0)<0)
		{
		F_v(0)=(-0.4);
		}
		else
		{
		F_v(0)=0;
		}
		
		//if(theta(1)<(-1.61))
		if(theta_dot(1)>0)
		{
		F_v(1)=0.2;
		}
		//else if(theta(1)>(-1.58))
		else if(theta_dot(1)<0)
		{
		F_v(1)=(-0.2);
		}
		else
		{
		F_v(1)=0;
		}
	tau = M * y + C * theta_dot + g_v + F_v;
//	tau = M * y + C * theta_dot + g_v;
//	std::cout<<"	g_v="<<g_v.transpose();
	tau(0)=-tau(0);//to youbot
}

brics_actuator::JointTorques generateJointTorqueMsg(Vector2d& joints)
{
  brics_actuator::JointTorques m_joint_torques;
  //Ros component negates torque values for joints with negative direction (all joints except joint 3)
//  joints(0) = -joints(0);
  std::stringstream jointName;
  m_joint_torques.torques.clear();
  Vector5d joint5;
  //joint5 << 0,0,joints(0),0,0;
  joint5 << 0,0,joints(0),joints(1),0;
//  joint5 << 0,0,joints,0;
 
 for (int i = 0; i < 5; i++)
  {
    brics_actuator::JointValue joint;
    joint.value = joint5(i);
    joint.unit = boost::units::to_string(boost::units::si::newton_meter);
    jointName.str("");
    jointName << "arm_joint_" << (i + 1);
    joint.joint_uri = jointName.str();

    m_joint_torques.torques.push_back(joint);
  }
  return m_joint_torques;
}

brics_actuator::JointPositions generateJointPositionMsg(Vector2d& joints)
{
  brics_actuator::JointPositions joint_position_msg;

  Vector5d joint5;
  joint5 << 0.05, 0.05, joints, 0.1107;


  std::stringstream jointName;
  joint_position_msg.positions.clear();

  for (int i = 0; i < 5; i++)
  {
    brics_actuator::JointValue joint;

    joint.value = joint5(i);
    joint.unit = boost::units::to_string(boost::units::si::radians);
    jointName.str("");
    jointName << "arm_joint_" << (i + 1);
    joint.joint_uri = jointName.str();

    joint_position_msg.positions.push_back(joint);
  }

  return joint_position_msg;
}

void writeToFile(std::ofstream & file, Vector2d & vect)
{
	for (int i = 0; i < vect.size(); i++)
	{
	file << vect[i] << "\t";
	}
	file << std::endl;
}

class Subscri
{
public:
	Vector2d actual_p;
	Vector2d actual_v;
	Subscri(Vector2d init_pv):actual_p(init_pv)
	{
		actual_v<<0,0;
	}

	~Subscri()
	{
	}

	void callback(const sensor_msgs::JointState& input)
	{
		for(int i=0;i<2;i++)
		{
			actual_p(i)=input.position[i+2];
			actual_v(i)=input.velocity[i+2];//!aufpassen, vorzeichen
		}
	}
};
void youbot2torque(Vector2d& theta_y,Vector2d& theta_y_dot,Vector2d& theta_c,Vector2d& theta_c_dot)
{
  theta_c(0) = joint_offsets[2] - theta_y(0);
  theta_c(1) = theta_y(1) - joint_offsets[3];
  theta_c_dot(0) = -theta_y_dot(0);
  theta_c_dot(1) = theta_y_dot(1);
}
void youbot2Impedanz(Vector2d& theta_y,Vector2d& theta_y_dot,Vector2d& theta_I,Vector2d& theta_I_dot)
{
  theta_I(0) = theta_y(0) - joint_offsets[2];
  theta_I(1) = theta_y(1) - joint_offsets[3];
  theta_I_dot(0) = theta_y_dot(0);
  theta_I_dot(1) = theta_y_dot(1);
} 
int safecheck(Vector2d &theta_I,Vector2d &theta_I_dot,Vector2d &X_actual)
{
		if(theta_I(0)<-0.6||theta_I(0)>2||(X_actual(1)<(-2.179)*X_actual(0)&&X_actual(1)>0)||(X_actual(1)<(-0.684)*X_actual(0)&&X_actual(1)<0)||theta_I(1)<-1.58||theta_I(1)>1.58)
		{
			if(theta_I(0)<-0.6)
			{
				ROS_ERROR("Grenze1:theta_I(0)=%f<-0.6",theta_I(0));					
			}
			else if(theta_I(0)>2)
			{
				ROS_ERROR("Grenze2:theta_I(0)=%f>2",theta_I(0));
			}
			else if((X_actual(1)<(-2.179)*X_actual(0)&&X_actual(1)>0))
			{
				ROS_ERROR("Grenze3:%f < (-2.179)*%f ",X_actual(1),X_actual(0));
			}
			else if((X_actual(1)<(-0.684)*X_actual(0)&&X_actual(1)<0))
			{
				ROS_ERROR("Grenze4: %f < (-0.684)*%f",X_actual(1),X_actual(0));
			}
			else if(theta_I(1)<-1.58)
			{
				ROS_ERROR("Grenze5:theta_I(1)= %f <-1.58",theta_I(1));	
			}
			else if(theta_I(1)>1.58)
			{
				ROS_ERROR("Grenze6:theta_I(1)= %f >1.58",theta_I(1));
			}
			std::cout<<"	ros_shutdown\n";
			return 0;
		}
		else if(theta_I(0)>1.8||(X_actual(1)<(-4.2686)*X_actual(0)&&X_actual(1)>0)||theta_I(1)>1.4)
		{
			if(theta_I_dot(0)>0.2||theta_I_dot(1)>0.2)
			{
				ROS_ERROR("Velocity error23");
				return 0;
			}
		}
		else if(theta_I(0)<-0.48||(X_actual(1)<(-0.5209)*X_actual(0)&&X_actual(1)<0)||theta_I(1)<-1.4)
		{
			if(theta_I_dot(0)<-0.2||theta_I_dot(1)<-0.2)
			{
				ROS_ERROR("Velocity error14");
				return 0;
			}
		}
		else
		{
			return 1;
		}
}

#endif
