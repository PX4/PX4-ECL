from sympy import *
from code_gen import *

def create_cov_matrix(i, j):
    if j >= i:
        return Symbol("P" + str(i) + str(j), real=True)
    else:
        return 0

def create_symmetric_cov_matrix():
    # define a symbolic covariance matrix
    P = Matrix(3,3,create_cov_matrix)

    for index in range(3):
        for j in range(3):
            if index > j:
                P[index,j] = P[j,index]

    return P

def create_symbol(name, real=True):
    symbol_name_list.append(name)
    return Symbol(name, real=True)

symbol_name_list = []

daz = create_symbol("daz", real=True) # IMU z axis delta angle measurement in body axes - rad
dvx = create_symbol("dvx", real=True) # IMU x axis delta velocity measurement in body axes - m/sec
dvy = create_symbol("dvy", real=True) # IMU y axis delta velocity measurement in body axes - m/sec
psi = create_symbol("psi", real=True)  # yaw angle of body frame wrt earth frame
vn = create_symbol("vn", real=True) # N velocity - m/sec
ve = create_symbol("ve", real=True) # E velocity - m/sec
dt = create_symbol("dt", real=True)  # dt (sec)
dazVar = create_symbol("dazVar", real=True) # IMU Z axis delta angle measurement variance (rad^2)
dvxVar = create_symbol("dvxVar", real=True) # IMU X axis delta velocity measurement variance (m/s)^2
dvyVar = create_symbol("dvyVar", real=True) # IMU Y axis delta velocity measurement variance (m/s)^2

# derive the body to nav direction transformation matrix
Tbn = Matrix([[cos(psi) , -sin(psi)],
              [sin(psi) ,  cos(psi)]])

# attitude update equation
psiNew = psi + daz

# velocity update equations
velNew = Matrix([vn,ve]) + Tbn*Matrix([dvx,dvy])

# Define the state vectors
stateVector = Matrix([vn,ve,psi])

# Define vector of process equations
newStateVector = Matrix([velNew,psiNew])

# Calculate state transition matrix
F = newStateVector.jacobian(stateVector)

# Derive the covariance prediction equations
# Error growth in the inertial solution is assumed to be driven by 'noise' in the delta angles and
# velocities, after bias effects have been removed.

# derive the control(disturbance) influence matrix from IMU noise to state noise
G = newStateVector.jacobian(Matrix([dvx,dvy,daz]))

# derive the state error matrix
distMatrix = Matrix([[dvxVar , 0 , 0],
                     [0 , dvyVar , 0],
                     [0 , 0 , dazVar]])

Q = G*distMatrix*G.T

# propagate covariance matrix
P = create_symmetric_cov_matrix()

P_new = F * P * F.T + Q

P_new_simple = cse(P_new, symbols("S0:1000"), optimizations='basic')

cov_prediction_code_generator = CodeGenerator("./covariance_prediction_generated.cpp")
cov_prediction_code_generator.print_string("Equations for covariance matrix prediction")
cov_prediction_code_generator.write_subexpressions(P_new_simple[0])
cov_prediction_code_generator.write_matrix(Matrix(P_new_simple[1]), "_ekf_gsf[model_index].P", True)
cov_prediction_code_generator.close()

# derive the covariance update equation for a NE velocity observation
velObsVar = create_symbol("velObsVar", real=True) # velocity observation variance (m/s)^2
H = Matrix([[1,0,0],
            [0,1,0]])

R = Matrix([[velObsVar , 0],
            [0 , velObsVar]])

S = H * P * H.T + R

S_simple = cse(S, symbols("S0:1000"), optimizations='basic')

innov_var_code_generator = CodeGenerator("./innovation_variance_generated.cpp")
innov_var_code_generator.print_string("Equations for NE velocity innovation variance")
innov_var_code_generator.write_subexpressions(S_simple[0])
innov_var_code_generator.write_matrix(Matrix(S_simple[1]), "S", True)
innov_var_code_generator.close()

# Calculate Kalman gain
K = (P * H.T) / S

K_simple = cse(K, symbols("SK0:1000"), optimizations='basic')

kalman_gain_code_generator = CodeGenerator("./Kalman_gain_generated.cpp")
kalman_gain_code_generator.print_string("Equations for NE velocity Kalman gain")
kalman_gain_code_generator.write_subexpressions(K_simple[0])
kalman_gain_code_generator.write_matrix(Matrix(K_simple[1]), "K", False)
kalman_gain_code_generator.close()

# Calculate updated covariance matrix
P_new = P - K*S*K.T

P_new_simple = cse(P_new, symbols("SP0:1000"), optimizations='basic')

cov_update_code_generator = CodeGenerator("./covariance_update_generated.cpp")
cov_update_code_generator.print_string("Equations for covariance matrix update")
cov_update_code_generator.write_subexpressions(P_new_simple[0])
cov_update_code_generator.write_matrix(Matrix(P_new_simple[1]), "_ekf_gsf[model_index].P", True)
cov_update_code_generator.close()