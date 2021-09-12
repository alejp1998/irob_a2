#! /usr/bin/env python3

"""
    # Alejandro Jarabo Peñas
    # aljp@kth.se
"""

import math

def scara_IK(point):
    x = point[0]
    y = point[1]
    z = point[2]

    # The length of the links is: 
    l = [0.07, 0.3, 0.35]
    

    # The forward kinematics are :
    # x = l_1c_1 + l_2c_12 + l_0 
    # y = l_1s_1 + l_2s_12
    # z = d_3

    # Therefore, solving for theta_1 and theta_2, we have that:
    c_2 = (x**2 + y**2 - (l[0]**2 + l[1]**2+ l[2]**2))/(2*l[1]*l[2])
    theta_2 = math.acos(c_2)
    s_2 = math.sin(theta_2)

    k_1 = l[1] + l[2]*c_2
    k_2 = l[2]*s_2

    theta_1 = math.atan2(y,x-l[0]) - math.atan2(k_2,k_1)

    # The position of the third joint coincides with the desired z value
    d_3 = z

    # So, the final positon config vector is:
    q = [theta_1, theta_2, d_3]

    return q

def kuka_IK(point, R, joint_positions):
    x = point[0]
    y = point[1]
    z = point[2]
    q = joint_positions #it must contain 7 elements

    """
    Fill in your IK solution here and return the seven joint values in q
    """

    return q
