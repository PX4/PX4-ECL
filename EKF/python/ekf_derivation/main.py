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
        return Symbol("P(" + str(i) + "," + str(j) + ")", real=True)
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

def generate_ccode(Input):
    f_output = open(Input["save_location"], "a+")
    
    write_string = ""
    
    if Input["shape"][0] > 0 and Input["shape"][1] > 0:
        # write a matrix
        write_string = "float " + Input["array_identifier"] + "[" + str(Input["shape"][0]) + "][" + str(Input["shape"][1]) + "] = {};\n"
        
        for j in range(0, Input["shape"][1]):
            for i in range(0, Input["shape"][0]):
                if j >= i or not Input["symetric_matrix"]:
                    write_string = write_string + Input["array_identifier"] + "(" + str(i) + "," + str(j) + ") = " + ccode(Input["data"][i,j]) + ";\n"
    elif  Input["shape"][0] == 0 and Input["shape"][1] == 0:
        for item in Input["data"]:
            write_string = write_string + "float " + str(item[0]) + " = " + ccode(item[1]) + ";\n"
    
    f_output.write(write_string)
    f_output.close()
            

def create_symbol(name, real=True):
    symbol_name_list.append(name)
    return Symbol(name, real=True)

def set_values_for_cov_update(Input):
    f_output = open(Input["save_location"], "a+")
    
    for index,item in enumerate(symbol_name_list):
        write_string = "float " + item + " = " + str(symbol_value_list[item]) + ";\n"
        f_output.write(write_string)
    
    f_output.close()
    
def set_values_for_cov_update_matlab(Input):
    f_output = open(Input["save_location"], "a+")
    
    for item in symbol_values_matlab.keys():
        write_string = "float " + item + " = " + str(symbol_values_matlab[item]) + ";\n"
        f_output.write(write_string)
    
    f_output.close()
    
def set_symbol_values():
    symbol_value_list = {
           "dt" : 0.01,
           "g"  : 9.81,
           "R_mag" : 0.01,
           "R_baro": 4.0,
           "R_hor_vel" : 0.1**2,
           "R_vert_vel" : 0.2**2,
           "R_hor_pos" : 0.5**2,
           "d_ang_x" : 0.1,
           "d_ang_y" : 0.1,
           "d_ang_z": 0.1,
           "d_v_x"  : 0.1,
           "d_v_y"  : 0.1,
           "d_v_z"  :0.2,
           "d_ang_x_var" : 0.01**2,
           "d_ang_y_var" : 0.01**2,
           "d_ang_z_var" : 0.01**2,
           "d_v_x_var"  : 0.1**2,
           "d_v_y_var"  : 0.1**2,
           "d_v_z_var"  : 0.1**2,
           "qw"     : 1,
           "qx"     : 0,
           "qy"     : 0,
           "qz"     : 0,
           "vx"     : 0.1,
           "vy"     : 0.2,
           "vz"     : 0.3,
           "px"     : 0,
           "py"     : 0,
           "pz"     : 0,
           "d_ang_bx" : 0.001,
           "d_ang_by" : 0.002,
           "d_ang_bz" : 0.003,
           "d_vel_bx" : 0.001,
           "d_vel_by" : 0.002,
           "d_vel_bz" : 0.003,
           "ix" : 0.2,
           "iy" : 0,
           "iz" : 0.5,
           "ibx" : 0.01,
           "iby" : 0.02,
           "ibz" : 0.03,
           "wx" : -1,
           "wy" : 2
            }
    return symbol_value_list

def set_symbol_values_matlab():
    symbol_value_list = {
           "dt" : 0.01,
           "g"  : 9.81,
           "R_mag" : 0.01,
           "R_baro": 4.0,
           "R_hor_vel" : 0.1**2,
           "R_vert_vel" : 0.2**2,
           "R_hor_pos" : 0.5**2,
           "dax" : 0.1,
           "day" : 0.1,
           "daz": 0.1,
           "dvx"  : 0.1,
           "dvy"  : 0.1,
           "dvz"  :0.2,
           "daxVar" : 0.01**2,
           "dayVar" : 0.01**2,
           "dazVar" : 0.01**2,
           "dvxVar"  : 0.1**2,
           "dvyVar"  : 0.1**2,
           "dvzVar"  : 0.1**2,
           "q0"     : 1,
           "q1"     : 0,
           "q2"     : 0,
           "q3"     : 0,
           "vx"     : 0.1,
           "vy"     : 0.2,
           "vz"     : 0.3,
           "px"     : 0,
           "py"     : 0,
           "pz"     : 0,
           "dax_b" : 0.001,
           "day_b" : 0.002,
           "daz_b" : 0.003,
           "dvx_b" : 0.001,
           "dvy_b" : 0.002,
           "dvz_b" : 0.003,
           "ix" : 0.2,
           "iy" : 0,
           "iz" : 0.5,
           "ibx" : 0.01,
           "iby" : 0.02,
           "ibz" : 0.03,
           "wx" : -1,
           "wy" : 2
            }
    return symbol_value_list
        

symbol_name_list = []
symbol_value_list = set_symbol_values()
symbol_values_matlab = set_symbol_values_matlab()

dt = create_symbol("dt", real=True)  # dt
g = create_symbol("g", real=True) # gravity constant

r_mag = create_symbol("R_mag", real=True)  # magnetometer measurement noise variance
r_baro = create_symbol("R_baro", real=True)    # barometer noise variance
r_hor_vel = create_symbol("R_hor_vel", real=True) # horizontal velocity noise variance
r_ver_vel = create_symbol("R_vert_vel", real=True) # vertical velocity noise variance
r_hor_pos = create_symbol("R_hor_pos", real=True) # horizontal position noise variance

# inputs, integrated gyro measurements
d_ang_x = create_symbol("dax", real=True)  # delta angle x
d_ang_y = create_symbol("day", real=True)  # delta angle y
d_ang_z = create_symbol("daz", real=True)  # delta angle z

d_ang = Matrix([d_ang_x, d_ang_y, d_ang_z])

# inputs, integrated accelerometer measurements
d_v_x = create_symbol("dvx", real=True)  # delta velocity x
d_v_y = create_symbol("dvy", real=True)  # delta velocity y
d_v_z = create_symbol("dvz", real=True)  # delta velocity z

d_v = Matrix([d_v_x, d_v_y,d_v_z])

u = Matrix([d_ang, d_v])

# input noise
d_ang_x_var = create_symbol("daxVar", real=True)
d_ang_y_var = create_symbol("dayVar", real=True)
d_ang_z_var = create_symbol("dazVar", real=True)

d_v_x_var = create_symbol("dvxVar", real=True)
d_v_y_var = create_symbol("dvyVar", real=True)
d_v_z_var = create_symbol("dvzVar", real=True)

var_u = Matrix.diag(d_ang_x_var, d_ang_y_var, d_ang_z_var, d_v_x_var, d_v_y_var, d_v_z_var)

# define state vector
    
# attitude quaternion
qw = create_symbol("q0", real=True)  # quaternion real part
qx = create_symbol("q1", real=True)  # quaternion x component
qy = create_symbol("q2", real=True)  # quaternion y component
qz = create_symbol("q3", real=True)  # quaternion z component

q = Matrix([qw,qx,qy,qz])
R_to_earth = quat2Rot(q)
R_to_body = R_to_earth.T


# velocity in NED local frame
vx = create_symbol("vx", real=True)  # north velocity
vy = create_symbol("vy", real=True)  # east velocity
vz = create_symbol("vz", real=True)  # down velocity

v = Matrix([vx,vy,vz])

# position in NED local frame
px = create_symbol("px", real=True)  # north position
py = create_symbol("py", real=True)  # east position
pz = create_symbol("pz", real=True)  # down position

p = Matrix([px,py,pz])

# delta angle bias
d_ang_bx = create_symbol("dax_b", real=True)  # delta angle bias x
d_ang_by = create_symbol("day_b", real=True)  # delta angle bias y
d_ang_bz = create_symbol("daz_b", real=True)  # delta angle bias z

d_ang_b = Matrix([d_ang_bx, d_ang_by, d_ang_bz])
d_ang_true = d_ang - d_ang_b


# delta velocity bias
d_vel_bx = create_symbol("dvx_b", real=True)  # delta velocity bias x
d_vel_by = create_symbol("dvy_b", real=True)  # delta velocity bias y
d_vel_bz = create_symbol("dvz_b", real=True)  # delta velocity bias z

d_vel_b = Matrix([d_vel_bx, d_vel_by, d_vel_bz])

d_vel_true = d_v - d_vel_b

# earth magnetic field vector
ix = create_symbol("ix", real=True)  # earth magnetic field x component
iy = create_symbol("iy", real=True)  # earth magnetic field y component
iz = create_symbol("iz", real=True)  # earth magnetic field z component

i = Matrix([ix,iy,iz])

# magnetometer bias in body frame
ibx = create_symbol("ibx", real=True)  # earth magnetic field bias in body x
iby = create_symbol("iby", real=True)  # earth magnetic field bias in body y
ibz = create_symbol("ibz", real=True)  # earth magnetic field bias in body z

ib = Matrix([ibx,iby,ibz])

# wind in local NE frame
wx = create_symbol("wx", real=True)  # wind in north direction
wy = create_symbol("wy", real=True)  # wind in east direction

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
Psimple = Matrix(P_new_simple[1])


code_gen_data = {
     "data" : P_new_simple[0],
     "shape" : (0,0),
     "save_location" : "./generated_python.cpp"
     }

generate_ccode(code_gen_data)

            
code_gen_data = {
     "data" : Psimple,
     "shape" : (24,24),
     "array_identifier" : "nextP",
     "symetric_matrix" : True,
     "save_location" : "./generated_python.cpp"
     }

generate_ccode(code_gen_data)

#set_values_for_cov_update({"save_location" : "./generated_python.cpp"})

#set_values_for_cov_update_matlab({"save_location" : "./generated_matlab.cpp"})


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

a= P_new.subs([(d_ang_x, 0.1), (d_ang_y, 0.2), (d_ang_z, 0.3), (d_ang_bx, 0.1), (d_ang_by, 0.2), (d_ang_bz, 0.3), (qw, 1), (qx, 0.3), (qy, 0.4), (qz, 0.5), (dt, 0.01), (d_v_x, 0.1), (d_v_y, 0.2), (d_v_z, 0.3), (d_vel_bx, 0.1), (d_vel_by, 0.2), (d_vel_bz, 0.3), (d_ang_x_var, 0.1), (d_ang_y_var, 0.2), (d_ang_z_var, 0.3), (d_v_x_var, 0.1), (d_v_y_var, 0.2), (d_v_z_var, 0.3)])


for index in range(24):
    for j in range(24):
        a = a.subs([(P[index,j], 0.1)])
       
        
sum_P = 0
for index in range(24):
    for j in range(24):
        sum_P = sum_P + a[index,j]
        
for index in range(24):
    sum_P = 0
    for j in range(24):
        sum_P = sum_P + a[index,j]
    
    print("sum %s %s") % (str(index), str(sum_P))
        
s = 0
for index in range(24):
        print(a[4,index])  
        
# counter number of operations
mult = 0
sumations = 0

for item in P_new_simple[0]:
    expression_string = str(item)
    
    for index,char in enumerate(expression_string):
        if char == "+" or char == "-":
            sumations = sumations + 1
            
        if char == "*":
            if index > 0 and expression_string[index-1] != "*":
                if index < len(expression_string) -1 and expression_string[index+1] != "*":
                    mult = mult +1
       
    
for item in P_new_simple[1][0]:
    expression_string = str(item)
    
    for index,char in enumerate(expression_string):
        if char == "+" or char == "-":
            sumations = sumations + 1
            
        if char == "*":
            if index > 0 and expression_string[index-1] != "*":
                if index < len(expression_string) -1 and expression_string[index+1] != "*":
                    mult = mult +1


f = open("./test.cpp", "r")

sum_matlab = 0
mult_matlab = 0

started = False

for line in f.readlines():
    if "begin" in line:
        started = True
        
    if "end" in line:
         break
        
    if started:
        for index,char in enumerate(line):
            if char == "+" or char == "-":
                sum_matlab = sum_matlab + 1
        
            if char == "*":
                if index > 0 and line[index-1] != "*":
                    if index < len(line) -1 and line[index+1] != "*":
                        mult_matlab = mult_matlab +1


f.close()

