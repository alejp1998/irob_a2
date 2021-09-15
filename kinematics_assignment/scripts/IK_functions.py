#! /usr/bin/env python3

"""
    # Alejandro Jarabo Peñas
    # aljp@kth.se
"""

from math import sin, asin, cos, acos, atan2, pi
import numpy as np
import tf2_ros as tf
import rospy

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
    c_2 = (x**2 - 2*l[0]*x + y**2 + l[0]**2 - l[1]**2 - l[2]**2)/(2*l[1]*l[2])
    theta_2 = acos(c_2)
    s_2 = sin(theta_2)

    k_1 = l[1] + l[2]*c_2
    k_2 = l[2]*s_2

    theta_1 = atan2(y,x-l[0]) - atan2(k_2,k_1)

    # The position of the third joint coincides with the desired z value
    d_3 = z

    # So, the final positon config vector is:
    q = [theta_1, theta_2, d_3]

    print('Desired cartesian pose: \n x = {}, y = {}, z = {}'.format(x,y,z))
    print('Desired joint pose: \n theta_1 = {}, theta_2 = {}, d_3 = {}\n'.format(theta_1,theta_2,d_3))

    return q

def kuka_IK(point, R, joint_positions):
    x = point[0]
    y = point[1]
    z = point[2]
    q = joint_positions #it must contain 7 elements

    # BASE LINK LENGTH
    B = 0.311
    # MIDDLE LINKS LENGTH
    L = 0.4
    M = 0.39
    # END EFFECTOR LENGTH
    E = 0.078

    # DH TABLE PARAMS
    alpha = [pi/2, -pi/2, -pi/2, pi/2, pi/2, -pi/2,   0]
    d     = [   0,     0,     L,    0,    M,     0,   0]
    a     = [   0,     0,     0,    0,    0,     0,   0]

    # Tolerance
    tol = 0.0001
    eps_x = np.ones(6)
    rate = rospy.Rate(100)

    #Desired x = [x, y, z, ...]
    R = np.array(R)
    theta = -asin(R[2,0])
    psi = atan2(R[2,1]/cos(theta), R[2,2]/cos(theta))
    phi = atan2(R[1,0]/cos(theta), R[0,0]/cos(theta))
    euler_desired = [psi, theta, phi]
    x_desired = np.hstack((np.array(point), euler_desired))

    while abs(max(eps_x)) > tol :
        # Create transformation matrices between two links 
        # according to Modified DH convention with given parameters  
        Ts = [np.array([
                [ cos(q[i]), -sin(q[i])*cos(alpha[i]),   sin(q[i])*sin(alpha[i]),      a[i]*cos(q[i])],
                [ sin(q[i]),  cos(q[i])*cos(alpha[i]),  -cos(q[i])*sin(alpha[i]),  a[i]*sin(alpha[i])],
                [         0,            sin(alpha[i]),             cos(alpha[i]),                d[i]],
                [         0,                        0,                         0,                   1]
         ]) for i in range(len(q)) ]
        
        #BASE LENGTH TRANSFORM 
        base_translation = np.eye(4)
        base_translation[2,3] = B
        # FINAL TRANSFORMS 
        TB1, T12, T23, T34, T45, T56, T67 = Ts[0], Ts[1], Ts[2], Ts[3], Ts[4], Ts[5], Ts[6]
        # END EFFECTOR TRANSFORM
        end_effector_translation = np.eye(4)
        end_effector_translation[2,3] = E

        # TRANSFORM FROM BASE TO END-EFFECTOR
        TB_st = np.dot(base_translation, TB1)
        TB2 = np.dot(TB_st, T12)
        TB3 = np.dot(TB2, T23)
        TB4 = np.dot(TB3, T34)
        TB5 = np.dot(TB4, T45)
        TB6 = np.dot(TB5, T56)
        TB7 = np.dot(TB6, T67)
        TBE = np.dot(TB7, end_effector_translation)

        # TRANSLATION VECTORS
        p = []
        p.append(TB_st[0:3, 3])
        p.append(TB2[0:3, 3])
        p.append(TB3[0:3, 3])
        p.append(TB4[0:3, 3])
        p.append(TB5[0:3, 3])
        p.append(TB6[0:3, 3])
        p.append(TB7[0:3, 3])
        p.append(TBE[0:3, 3])

        # Z ROTATION VECTORS
        z = []
        z.append(TB_st[0:3, 2])
        z.append(TB2[0:3, 2])
        z.append(TB3[0:3, 2])
        z.append(TB4[0:3, 2])
        z.append(TB5[0:3, 2])
        z.append(TB6[0:3, 2])
        z.append(TB7[0:3, 2])
        z.append(TBE[0:3, 2])

        # COMPUTE EULER ANGLES: psi, theta and phi
        # Rotation matrix
        R = TBE[0:3, 0:3]

        # Euler angles estimation
        theta = -asin(R[2,0])
        psi = atan2(R[2,1]/cos(theta), R[2,2]/cos(theta))
        phi = atan2(R[1,0]/cos(theta), R[0,0]/cos(theta))

        euler_hat = [psi, theta, phi]

        # Position estimation
        x_hat = np.hstack((np.array([TBE[0][3], TBE[1][3], TBE[2][3]]),euler_hat))
        # Calculate position error
        eps_x = x_hat - x_desired

        # Calculate Jacobian
        J = []
        for i in range(len(q)) : 
            J_temp = np.cross(z[i], (p[-1] - p[i]))
            J.append(np.hstack((J_temp, z[i])))
        J = np.array(J)

        # Calculate pose error
        J_pseudo_inv = np.linalg.pinv(J.T)
        eps_theta = np.matmul(J_pseudo_inv, eps_x)

        # Obtain new pose estimation
        q = q - eps_theta
        print('New q: ')
        print(q)
        rate.sleep()


    return q
