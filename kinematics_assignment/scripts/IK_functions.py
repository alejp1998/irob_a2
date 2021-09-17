#! /usr/bin/env python3

"""
    # Alejandro Jarabo
    # aljp@kth.se
"""

from math import sin, asin, cos, acos, atan2, pi
import numpy as np

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
    tol = 0.001

    #Desired pose from R and point
    R = np.array(R)
    n_d = np.array(R[:3,0])
    o_d = np.array(R[:3,1])
    a_d = np.array(R[:3,2])
    x_d = np.array(point)
    
    not_small_enough_eps = True
    while not_small_enough_eps :
        not_small_enough_eps = False

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
        # END EFFECTOR TRANSFORM
        end_effector_translation = np.eye(4)
        end_effector_translation[2,3] = E

        # FINAL TRANSFORM MATRIX AND JACOBIAN PARAMETERS
        # Jacobian Params
        pose_p = np.zeros(shape = (7, 3, 1))
        J_o = np.zeros(shape = (7, 3, 1))
        J_p = np.zeros(shape = (7, 3, 1))
        J = np.zeros(shape = (6, 7))
        # Transform from 0 to 7
        T07 = np.eye(4) 
        for i in range(len(q)) :
            J_o[i] = T07[:3, [2]]
            pose_p[i] = T07[:3, [3]]
            T07 = np.dot(T07,Ts[i])
        # Transform from 0 to End-Effector
        T0E = np.dot(T07,end_effector_translation)
        # Transform from Base to End-Effector
        TBE = np.dot(base_translation,T0E)
        
        # End-Effector Transformed Rotation
        n_e = TBE[:3, [0]]
        o_e = TBE[:3, [1]]
        a_e = TBE[:3, [2]]
        
        # Computation of Jacobian and Pseudo-Inverse Jacobian
        pose_e = T0E[:3, [3]]
        for i in range(len(q)) :
            J_p[i] = np.transpose(np.cross(np.transpose(J_o[i]), np.transpose(pose_e - pose_p[i])))
        for i in range(len(q)) :
            for j in range(len(point)) :
                J[j][i] = J_p[i][j][0]
                J[j+3][i] = J_o[i][j][0]
        
        J_pseudo_inv = np.dot(np.transpose(J), (np.linalg.inv(np.dot(J, np.transpose(J)))))

        # Compute current pose
        x_c = np.array(TBE[:3, 3])

        # Position error
        eps_x = x_c - x_d
        eps_x = np.array([[eps_x[0]], [eps_x[1]],[eps_x[2]]])

        # Orientation error
        eps_n = np.cross(np.transpose(n_e), np.transpose(n_d))
        eps_o = np.cross(np.transpose(o_e), np.transpose(o_d))
        eps_a = np.cross(np.transpose(a_e), np.transpose(a_d))
        eps_o = 0.5*np.transpose(eps_n + eps_o + eps_a)

        # Complete pose error
        eps_pose = np.append(eps_x, eps_o, axis=0)

        # Joint angle error from pose error
        eps_joint_angle = np.dot(J_pseudo_inv, eps_pose)

        # New q values
        for i in range(len(q)) : 
            q[i] = float(simplify_angle(q[i] - eps_joint_angle[i][0]))

        # Check if eps is smaller than tolerance
        pose_error = np.linalg.norm(eps_pose)
        print('Pose error: {}'.format(pose_error))
        if pose_error > tol : 
            not_small_enough_eps = True

        print('New q: ')
        print(q)

    print('Final q : ')
    print(q)
    
    return q

def simplify_angle(theta) :
    while (theta >= pi) or (theta < -pi) :
        if theta >= pi : 
            theta = theta - 2*pi
        elif theta < -pi :
            theta = theta + 2*pi

    return theta


# Test IK function
kuka_IK([-0.123, 0., 0.42], [[0, 0, -1], [0, 1, 0], [1, 0, 0]], [0, 1.11, 0, 1.1, 0, 1.52, 0])

