
#QFT Quatnum addition Algorithm Syed Affan Hussain 
import numpy as np
from numpy import pi
from qiskit import *


#_____________________________________________________________________________________________________________________________________________________________________________________________________________________                                                                    
def QFT_a(Qcirc,register_a,n):           #|quantumcirc.cp is used as Controlled R(k)gate or Phase rotation
                                         #|Qcirc.cp(theta,label ,control_state) is the Conrolled phase 
 Qcirc.h(register_a[n])                  #|rotation gate(the matrix representation is same as the Controlled R(k)).The input parameters the 
 for i in range(0, n):                   #|theta angle or phase angle,control qubit parameter and targeted qubit parameter. 
        theta = pi/float(2**(i+1))       #|What this function does is that it initiate a sequence  of hadamard gates and control gates
        c_qubit = register_a[n-(i+1)]    #|with decrease in order of control gates with each iteration decrease which means that 
        t_qubit =  register_a[n]         #|if the function with loops is used it will place a Hadamard gate on qubit n then it will perform a sequence of 
        Qcirc.cp(theta,t_qubit ,c_qubit) #|tasks afterwards n = 6(but it will be a dynamic variable(you will understand that in the later code)) the control qubit
                                         #|will be at 5 and it would change the its place of qubit every time the function replenishes ,the targeted qubit 
                                         #|will at n-1 qubit because it will change with each iteration(so in other terms with each iteration the target qubits will be set with decreasing order)
                                         #|This whole function is applied on register a which holds the  binary string of first integer as Qubits.
#__________________________________________________________________________________________________________________________________________________________________________________________________________________
def QFT_aplusb(Qcirc,register_a,register_b,n,factor): #|What this fucntion does is that is applies 
    l = len(register_b)
    
    for i in range(0, n+1):  
        if (n - i) > l - 1:   
            pass
        else:                                       #|
         theta = (factor*pi)/float(2**(i))                  #|
         c_qubit = register_a[n]                    #|
         t_qubit =  register_b[n-i]                 #|
         Qcirc.cp(theta, t_qubit, c_qubit)          #|
#__________________________________________________________________________________________________________________________________________________________________________________________________________________
#_____________________________________________________________________________________________________________________________________________________________________________________________________________________                                                                    
def QFT_inverse_a(Qcirc,register_a,n):        #|quantumcirc.cp is used as Controlled R(k)gate or Phase rotation
                                         #|Qcirc.cp(theta,label ,control_state) is the Conrolled phase 
                                         #|rotation gate(the matrix representation is same as the Controlled R(k)).The input parameters the 
    for i in range(0, n):                   #|theta angle or phase angle,control qubit parameter and targeted qubit parameter. 
        theta = (-1*pi)/float(2**(n-1))       #|What this function does is that it initiate a sequence  of hadamard gates and control gates
        c_qubit = register_a[n]    #|with decrease in order of control gates with each iteration decrease which means that 
        t_qubit =  register_a[i]         #|if the function with loops is used it will place a Hadamard gate on qubit n then it will perform a sequence of 
        Qcirc.cp(theta,t_qubit ,c_qubit) #|tasks afterwards n = 6(but it will be a dynamic variable(you will understand that in the later code)) the control qubit
    Qcirc.h(register_a[n])                  #|will be at 5 and it would change the its place of qubit every time the function replenishes ,the targeted qubit 
                                         #|will at n-1 qubit because it will change with each iteration(so in other terms with each iteration the target qubits will be set with decreasing order)
                                        #|This whole function is applied on register a which holds the  binary string of first integer as Qubits.
#__________________________________________________________________________________________________________________________________________________________________________________________________________________

def Quantum_add(reg_a, reg_b, Qcirc,factor ):
    
    #This function performs the QFT Adder circuit
    n = len(reg_a) - 1
    #Compute the Fourier transform of register a
    for i in range(0, n+1):
        QFT_a(Qcirc, reg_a, n-i)
    Qcirc.barrier()
    #change state of 2nd register w.r.t register
    for i in range(0, n+1):
        QFT_aplusb(Qcirc, reg_a , reg_b, n-i,factor)
    #Compute the inverse Fourier transform of register a
    Qcirc.barrier()
    for i in range(0, n + 1):
       QFT_inverse_a(Qcirc, reg_a, i)
    Qcirc.barrier()
    
 

#__________________________________________________________________________________________________________________________________________________________________________________________________________________
#__________________________________________________________________________________________________________________________________________________________________________________________________________________
def bin_to_dec(b):
    return int(b, 2)
 
def dec_to_bin(u):
    r = bin(u)[2:]
    return r
#_____________________________________________________________________________________________________________________________________________________________________________________________________________________ 
bs1 = dec_to_bin(int(input("Enter First Integer (Less than 500)to add: "))) #|This is done because Qiskit simulation only perform simulation with 32 qubit
bs2 = dec_to_bin(int(input("Enter Second Integer(Less than 500) to add: ")))#|and keep in mind that (n +1 ) of qubits are required for first and second binary string 
#bs means binary string                                                     #|and for classical storage register in orde to store output so 9 qubits should be encouraged to perform operation
lenbs1 = len(bs1)
lenbs2 = len(bs2)
                                                             #|1.Here length of both Binary strings are 
if (lenbs1 or lenbs2 )> 9 :                                  #|compared beacause Quatntum fourier transform of same length of binary string 
    print("The integer you input is greater that 9 bits ") #|helps better in computation rather than differnt length binary strings.
else:                                                      #|2.The first IF checks if the lenth of string if the length is greater  
  if lenbs1 != lenbs2:                                     #|than 9 it will throw an error
    if lenbs1 > lenbs2:
         bs1, bs2= bs2, bs1
    lenbs2, lenbs1 = lenbs1, lenbs2
reg_a = QuantumRegister(lenbs1,"reg_a") #Holds the first number
reg_b = QuantumRegister(lenbs2,"reg_b") #Holds the second number
ancila = QuantumRegister(lenbs1+lenbs2,"ancila")
classic_reg = ClassicalRegister(lenbs1+lenbs2,"classic_reg") #Holds the final output
x = QuantumRegister(1)
Qcirc = QuantumCircuit(ancila,reg_b, reg_a, x,classic_reg,name="Qcirc")#
Qcirc.x(x)
for i in range(0, lenbs1):
        if bs1[i] == "1":
            Qcirc.x(reg_a[lenbs1-(i+1)])
    #Flip the corresponding qubit in register b if a bit in the 
    #string second is a 1
for i in range(0, lenbs1):
        if bs2[i] == "1":
            Qcirc.x(reg_b[lenbs2-(i+1)])
mx = '1'
while(int(mx)!= 0):
    Quantum_add(ancila, reg_a, Qcirc, 1)
    Quantum_add(reg_b, x, Qcirc, -1)
    for i in range(len(reg_b)):
        Qcirc.measure(reg_b[i], classic_reg[i])
    result = execute(Qcirc, backend=Aer.get_backend('qasm_simulator'),
                    shots=2).result().get_counts(Qcirc.name)
    mx = list(result.keys())[0]
Qcirc.measure( classic_reg)
result = execute(Qcirc, backend=Aer.get_backend('qasm_simulator'),
            shots=2).result().get_counts(Qcirc.name)
#counts = counts = result.get_counts()
print(Qcirc.draw())
ans = int(next(iter(result)), 2)
print(ans)







