from sympy import *

# q: quaternion describing rotation from frame 1 to frame 2
# returns a rotation matrix derived form q which describes the same
# rotation
def quat2Rot(q):
    q0 = q[0]
    q1 = q[1]
    q2 = q[2]
    q3 = q[3]

    Rot = Matrix([[q0**2 + q1**2 - q2**2 - q3**2, 2*(q1*q2 - q0*q3), 2*(q1*q3 + q0*q2)],
                  [2*(q1*q2 + q0*q3), q0**2 - q1**2 + q2**2 - q3**2, 2*(q2*q3 - q0*q1)],
                   [2*(q1*q3-q0*q2), 2*(q2*q3 + q0*q1), q0**2 - q1**2 - q2**2 + q3**2]])

    return Rot

def create_cov_matrix(i, j):
    if j >= i:
        return Symbol("P_" + str(i) + "_" + str(j) + "", real=True)
    else:
        return 0
    
def quat_mult(p,q):
    r = Matrix([p[0] * q[0] - p[1] * q[1] - p[2] * q[2] - p[3] * q[3],
                p[0] * q[1] + p[1] * q[0] + p[2] * q[3] - p[3] * q[2],
                p[0] * q[2] - p[1] * q[3] + p[2] * q[0] + p[3] * q[1],
                p[0] * q[3] + p[1] * q[2] - p[2] * q[1] + p[3] * q[0]])
    
    return r

def create_symmetric_cov_matrix():
    # define a symbolic covariance matrix
    P = Matrix(24,24,create_cov_matrix)

    for index in range(24):
        for j in range(24):
            if index > j:
                P[index,j] = P[j,index]

    return P

dt = Symbol("dt", real=True)  # dt
g = Symbol("g", real=True) # gravity constant

r_mag = Symbol("R_mag", real=True)  # magnetometer measurement noise variance
r_baro = Symbol("R_baro", real=True)    # barometer noise variance
r_hor_vel = Symbol("R_hor_vel", real=True) # horizontal velocity noise variance
r_ver_vel = Symbol("R_vert_vel", real=True) # vertical velocity noise variance
r_hor_pos = Symbol("R_hor_pos", real=True) # horizontal position noise variance

# inputs, integrated gyro measurements
d_ang_x = Symbol("d_ang_x", real=True)  # delta angle x
d_ang_y = Symbol("d_ang_y", real=True)  # delta angle y
d_ang_z = Symbol("d_ang_z", real=True)  # delta angle z

d_ang = Matrix([d_ang_x, d_ang_y, d_ang_z])

# inputs, integrated accelerometer measurements
d_v_x = Symbol("d_v_x", real=True)  # delta velocity x
d_v_y = Symbol("d_v_y", real=True)  # delta velocity y
d_v_z = Symbol("d_v_z", real=True)  # delta velocity z

d_v = Matrix([d_v_x, d_v_y,d_v_z])

u = Matrix([d_ang, d_v])

# input noise
d_ang_x_var = Symbol("d_ang_x_var", real=True)
d_ang_y_var = Symbol("d_ang_y_var", real=True)
d_ang_z_var = Symbol("d_ang_z_var", real=True)

d_v_x_var = Symbol("d_v_x_var", real=True)
d_v_y_var = Symbol("d_v_y_var", real=True)
d_v_z_var = Symbol("d_v_z_var", real=True)

var_u = Matrix.diag(d_ang_x_var, d_ang_y_var, d_ang_z_var, d_v_x_var, d_v_y_var, d_v_z_var)

# define state vector
    
# attitude quaternion
qw = Symbol("qw", real=True)  # quaternion real part
qx = Symbol("qx", real=True)  # quaternion x component
qy = Symbol("qy", real=True)  # quaternion y component
qz = Symbol("qz", real=True)  # quaternion z component

q = Matrix([qw,qx,qy,qz])
R_to_earth = quat2Rot(q)
R_to_body = R_to_earth.T


# velocity in NED local frame
vx = Symbol("vx", real=True)  # north velocity
vy = Symbol("vy", real=True)  # east velocity
vz = Symbol("vz", real=True)  # down velocity

v = Matrix([vx,vy,vz])

# position in NED local frame
px = Symbol("px", real=True)  # north position
py = Symbol("py", real=True)  # east position
pz = Symbol("pz", real=True)  # down position

p = Matrix([px,py,pz])

# delta angle bias
d_ang_bx = Symbol("d_ang_bx", real=True)  # delta angle bias x
d_ang_by = Symbol("d_ang_by", real=True)  # delta angle bias y
d_ang_bz = Symbol("d_ang_bz", real=True)  # delta angle bias z

d_ang_b = Matrix([d_ang_bx, d_ang_by, d_ang_bz])
d_ang_true = d_ang - d_ang_b


# delta velocity bias
d_vel_bx = Symbol("d_vel_bx", real=True)  # delta velocity bias x
d_vel_by = Symbol("d_vel_by", real=True)  # delta velocity bias y
d_vel_bz = Symbol("d_vel_bz", real=True)  # delta velocity bias z

d_vel_b = Matrix([d_vel_bx, d_vel_by, d_vel_bz])

d_vel_true = d_v - d_vel_b

# earth magnetic field vector
ix = Symbol("ix", real=True)  # earth magnetic field x component
iy = Symbol("iy", real=True)  # earth magnetic field y component
iz = Symbol("iz", real=True)  # earth magnetic field z component

i = Matrix([ix,iy,iz])

# magnetometer bias in body frame
ibx = Symbol("ibx", real=True)  # earth magnetic field bias in body x
iby = Symbol("iby", real=True)  # earth magnetic field bias in body y
ibz = Symbol("ibz", real=True)  # earth magnetic field bias in body z

ib = Matrix([ibx,iby,ibz])

# wind in local NE frame
wx = Symbol("wx", real=True)  # wind in north direction
wy = Symbol("wy", real=True)  # wind in east direction

w = Matrix([wx,wy])

# state vector at arbitrary time t
state = Matrix([q,v,p,d_ang_b,d_vel_b,i,ib,w])

# define state propagation
q_new = quat_mult(q, Matrix([1, 0.5 * d_ang_true[0],  0.5 * d_ang_true[1],  0.5 * d_ang_true[2]]))

v_new = v + R_to_earth * d_vel_true + Matrix([0,0,g]) * dt

p_new = p + v * dt

d_ang_b_new = d_ang_b
d_vel_b_new = d_vel_b
i_new = i
ib_new = ib
w_new = w

# predicted state vector at time t + dt
state_new = Matrix([q_new, v_new, p_new, d_ang_b_new, d_vel_b_new, i_new, ib_new, w_new])

# state transition matrix
A = state_new.jacobian(state)

# B
G = state_new.jacobian(u)

P = create_symmetric_cov_matrix()

# propagate covariance matrix
P_new = A * P * A.T + G * var_u * G.T

for index in range(24):
    for j in range(24):
        if index > j:
            P_new[index,j] = 0


P_new_simple = cse(P_new, symbols("PS0:400"), optimizations='basic')


# magnetometer fusion
m_mag = R_to_body * i + ib

H_x_mag = Matrix([m_mag[0]]).jacobian(state)
H_y_mag = Matrix([m_mag[1]]).jacobian(state)
H_z_mag = Matrix([m_mag[2]]).jacobian(state)

K_x_mag = P * H_x_mag.T / (H_x_mag * P * H_x_mag.T + Matrix([r_mag]))
K_y_mag = P * H_y_mag.T / (H_y_mag * P * H_y_mag.T + Matrix([r_mag]))
K_z_mag = P * H_z_mag.T / (H_z_mag * P * H_z_mag.T + Matrix([r_mag]))

K_simple_x = cse(K_x_mag, symbols('KS0:200'))
K_simple_y = cse(K_y_mag, symbols('KS0:200'))
K_simple_z = cse(K_z_mag, symbols('KS0:200'))

# velocity fusion
m_v = v
H_v = m_v.jacobian(state)
K_v_x = P * H_v[0,:].T / (H_v[0,:] * P * H_v[0,:].T + Matrix([r_hor_vel]))
K_v_y = P * H_v[1,:].T / (H_v[1,:] * P * H_v[1,:].T + Matrix([r_hor_vel]))
K_v_z = P * H_v[2,:].T / (H_v[2,:] * P * H_v[2,:].T + Matrix([r_ver_vel]))

# position fusion
m_p = p
H_p = m_p.jacobian(state)

K_p_x = P * H_p[0,:].T / (H_p[0,:] * P * H_p[0,:].T + Matrix([r_hor_pos]))
K_p_y = P * H_p[1,:].T / (H_p[1,:] * P * H_p[1,:].T + Matrix([r_hor_pos]))
K_p_z = P * H_p[2,:].T / (H_p[2,:] * P * H_p[2,:].T + Matrix([r_baro]))

#a= P_new.subs([(d_ang_x, 0.1), (d_ang_y, 0.2), (d_ang_z, 0.3), (d_ang_bx, 0.1), (d_ang_by, 0.2), (d_ang_bz, 0.3), (qw, 1), (qx, 0.3), (qy, 0.4), (qz, 0.5), (dt, 0.01), (d_v_x, 0.1), (d_v_y, 0.2), (d_v_z, 0.3), (d_vel_bx, 0.1), (d_vel_by, 0.2), (d_vel_bz, 0.3), (d_ang_x_var, 0.1), (d_ang_y_var, 0.2), (d_ang_z_var, 0.3), (d_v_x_var, 0.1), (d_v_y_var, 0.2), (d_v_z_var, 0.3)])
#
#
#for index in range(24):
#    for j in range(24):
#        a = a.subs([(P[index,j], 0.1)])
#       
#        
#sum_P = 0
#for index in range(24):
#    for j in range(24):
#        sum_P = sum_P + a[index,j]
#        
#for index in range(24):
#    sum_P = 0
#    for j in range(24):
#        sum_P = sum_P + a[index,j]
#    
#    print("sum %s %s") % (str(index), str(sum_P))
#        
#s = 0
#for index in range(24):
#        print(a[4,index])  
#        
## counter number of operations
#mult = 0
#sumations = 0
#
#for item in P_new_simple[0]:
#    expression_string = str(item)
#    
#    for index,char in enumerate(expression_string):
#        if char == "+" or char == "-":
#            sumations = sumations + 1
#            
#        if char == "*":
#            if index > 0 and expression_string[index-1] != "*":
#                if index < len(expression_string) -1 and expression_string[index+1] != "*":
#                    mult = mult +1
#       
#    
#for item in P_new_simple[1][0]:
#    expression_string = str(item)
#    
#    for index,char in enumerate(expression_string):
#        if char == "+" or char == "-":
#            sumations = sumations + 1
#            
#        if char == "*":
#            if index > 0 and expression_string[index-1] != "*":
#                if index < len(expression_string) -1 and expression_string[index+1] != "*":
#                    mult = mult +1
#
#
#f = open("./test.cpp", "r")
#
#sum_matlab = 0
#mult_matlab = 0
#
#started = False
#
#for line in f.readlines():
#    if "begin" in line:
#        started = True
#        
#    if "end" in line:
#         break
#        
#    if started:
#        for index,char in enumerate(line):
#            if char == "+" or char == "-":
#                sum_matlab = sum_matlab + 1
#        
#            if char == "*":
#                if index > 0 and line[index-1] != "*":
#                    if index < len(line) -1 and line[index+1] != "*":
#                        mult_matlab = mult_matlab +1
#
#
#f.close()

