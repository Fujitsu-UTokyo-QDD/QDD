# -*- coding: utf-8 -*-
from qiskit import QuantumCircuit as QiskitCircuit
import sympy.ntheory.factor_ , numpy
from cmath import exp
import math
import collections

class QuantumCircuitForFactoringWithQulacs:
    """
    Shorのアルゴリズムの周期発見アルゴリズムのうちfの作用を汎用回路を作成する。
    量子ビットとして以下を用意
    ・第1量子ビット列:  mビット(m=2n)
    ・第2量子ビット列:  nビット
    ・補助ビットR1:     1ビット
    ・補助ビットR2:     nビット
    ・補助ビットR3:     1ビット
    ・繰り上がりビット: n-1ビット
    """
    def __init__(self, N, add_construct="Q-ADD", decompose_toffoli=False, use_clean_manegement=True):
        self.N = N
        self.a = 2
        self.n_bits = len(bin(N)) - 2
        self.m = 2*self.n_bits
        self.add_construct = add_construct
        self.decompose_toffoli=decompose_toffoli
        self.use_clean_management=use_clean_manegement
        self.index_i = 0

        self.num_qubits = 1
        if add_construct == "R-ADD":
            self.num_qubits = 5*self.n_bits + 1
        elif add_construct == "GT-ADD" or add_construct == "Q-ADD" or add_construct == "MIX-ADD":
            self.num_qubits = 4*self.n_bits + 2
        else:
            print("Error in add_construct")

        self.num_c2_not_gate = 0
        self.num_c_not_gate = 0
        self.num_not_gate = 0
        self.num_h_gate = 0
        self.num_t_gate = 0
        self.num_tdag_gate = 0
        self.num_C2R_gate = 0
        self.num_CR_gate = 0
        self.num_R_gate = 0
        self.num_single_gates = 0
        self.num_elementary_gates = 0
        self.qisc = QiskitCircuit(self.num_qubits)
        self.x = [i for i in range(self.m)]
        self.y = [i + len(self.x) for i in range(self.n_bits)]
        self.R3 = len(self.x) + len(self.y)
        self.R2 = [i + len(self.x) + len(self.y) + 1 for i in range(self.n_bits)]
        self.R1 = len(self.x) + len(self.y) + 1 + len(self.R2)

        if add_construct == "R-ADD":
            self.carry = [i + len(self.x) + len(self.y) + 1 + len(self.R2) + 1 for i in range(self.n_bits-1)]

        self.clean_ancilla = [i for i in range(self.num_qubits)]
        self.dirty_ancilla = []

    def compute_number_of_gates(self):
        self.num_elementary_gates = self.num_c2_not_gate+self.num_c_not_gate+self.num_not_gate+self.num_h_gate+self.num_t_gate+self.num_tdag_gate+14*self.num_C2R_gate+6*self.num_CR_gate+self.num_R_gate+self.num_single_gates

    def print_num_gates(self):
        self.compute_number_of_gates()
        print("# of C2-NOT:\t", self.num_c2_not_gate)
        print("# of C-NOT:\t", self.num_c_not_gate)
        print("# of not:\t", self.num_not_gate)
        print("# of H:\t\t", self.num_h_gate)
        print("# of T\t\t", self.num_t_gate)
        print("# of Tdg\t", self.num_tdag_gate)
        print("# of C2R gate\t", self.num_C2R_gate)
        print("# of CR gate\t", self.num_CR_gate)
        print("# of R gate\t", self.num_R_gate)
        print("--")
        print("# of gates\t", self.num_elementary_gates)
        print("# depth\t\t", self.circuit_depth)
        print("# of qubits\t", self.num_qubits)

    def to_dirty(self, number):
        if number in self.clean_ancilla:
            self.clean_ancilla.remove(number)
            self.dirty_ancilla.append(number)
    
    def to_clean(self, number):
        if number in self.dirty_ancilla:
            self.dirty_ancilla.remove(number)
            self.clean_ancilla.append(number)

    def find_order_algorithm(self, a, number=10000, measure=True):
        
        self.a = a
        self.add_X_gate(self.y[0])
        
        if self.add_construct == "Q-ADD":
            self.add_qft_circuit(self.n_bits+1, self.R2[0], False)
            
        self.add_mod_exp_circuit(a, True)

        if self.add_construct == "Q-ADD" or self.add_construct == "MIX-ADD":
            self.add_inverse_qft_circuit(self.n_bits+1, self.R2[0], False)
            
        self.add_inverse_qft_circuit(self.m)

        if measure:
            return self.measure_and_compute_order(number)

    def measure_and_compute_order(self, number):
        self.print_num_gates()
        state = QuantumState(self.num_qubits)

        ls = []
        for i in sampled_values:
            ls.append(int(format(i, 'b')[-self.m:], 2))

        multiple = 2**self.m

        # orderの計算
        counter = collections.Counter(ls)
        threshold = number * 0.005
        for key in counter.keys():
            if key !=0 and int(counter[key]) > threshold:
                if int(key) < multiple:
                    multiple = int(key)
        order = round(2**self.m / multiple)

        print()
        print(counter)
        print("2^m/r =", multiple, ", r =", order)

        return order

    def compute_prime_factor(self, order):
        tmp = 1
        for i in range(order//2):
            tmp = tmp*self.a % self.N
            p = math.gcd(tmp+1, self.N)
            q = int(self.N/p)
            factors = sorted([p, q])
        return factors

    def add_mod_exp_circuit(self, a, with_Hadamard=False):
        phi = sympy.ntheory.factor_.totient(self.N)
        for i in range(len(self.x)):
            self.index_i = i # MIX-ADDのため計算位置を記憶

            # 2**i を計算(a^phi = 1 mod Nを使う)
            tmp = 1
            for j in range(i):
                tmp = (tmp * 2) % phi

            # a**(2**i) mod N を計算
            tmp2 = 1
            for j in range(tmp):
                tmp2 = (tmp2 * a) % self.N 

            if with_Hadamard:
                self.add_H_gate(self.x[i])

            self.add_mod_mul_circuit(tmp2, [self.x[i]])

    def add_H_gate(self, target_bit):
        #self.circuit.add_H_gate(target_bit)
        self.qisc.h(target_bit)
        self.num_h_gate += 1
        self.to_dirty(target_bit)

    def add_mod_mul_circuit(self, d, controlled_bit_list=[]):

        # MIX-ADDでR-ADDからQ-ADDに切り替わる位置ではQFTを実施
        if self.add_construct == "MIX-ADD" and self.index_i == self.m - self.n_bits + 1:
            self.add_qft_circuit(self.n_bits+1, self.R2[0], False)

        self.add_mod_ps_circuit(d, controlled_bit_list, self.clean_ancilla, self.dirty_ancilla)
        if self.add_construct == "Q-ADD" or (self.add_construct == "MIX-ADD" and self.index_i >= self.m - self.n_bits + 1):
            self.add_inverse_qft_circuit(self.n_bits+1, self.R2[0], False)
        for i in range(self.n_bits):
            self.add_controlled_swap(controlled_bit_list, self.y[i], self.R2[i])
        if self.add_construct == "Q-ADD" or (self.add_construct == "MIX-ADD" and self.index_i >= self.m - self.n_bits + 1):
            self.add_qft_circuit(self.n_bits+1, self.R2[0], False)
        self.add_mod_ps_circuit((-pow(d, -1, self.N)) % self.N, controlled_bit_list, self.clean_ancilla, self.dirty_ancilla)

        for i in self.R2:
            self.to_clean(i)

    def add_controlled_swap(self, controlled_bit_list, qubit1, qubit2):

        # clean/dirtyはswapできず、swap対象のどちらのビットもdirtyになる
        ## controlledビットによってswapされない可能性もあるため

        # c-swap
        tmp_controlled_bit_list = [i for i in controlled_bit_list]
        if len(controlled_bit_list) > 1:
            print("error in add_controlled_swap: The number of controlled bits is over 2")
        self.add_ck_not_circuit([qubit2], qubit1)
        tmp_controlled_bit_list.append(qubit1)
        self.add_ck_not_circuit(tmp_controlled_bit_list, qubit2)
        self.add_ck_not_circuit([qubit2], qubit1)

    def add_mod_ps_circuit(self, d, controlled_bit_list=[], clean_bit_list=[], dirty_bit_list=[]):
        tmp_controlled_bit_list = [i for i in controlled_bit_list]
        for i in range(len(self.y)):
            tmp_controlled_bit_list.append(self.y[i])
            tmp = 1
            for j in range(i):
                tmp *= 2 % self.N
            tmp = tmp*d % self.N
            self.add_mod_add_circuit(tmp, tmp_controlled_bit_list)

            tmp_controlled_bit_list.remove(self.y[i])

    def add_mod_add_circuit(self, d, controlled_bit_list=[]):
        if self.add_construct == "GT-ADD":
            self.add_mod_add_circuit_with_GT_ADD(d, controlled_bit_list)
        else:
            assert(0)

    def add_mod_add_circuit_with_GT_ADD(self, d, controlled_bit_list=[]):
        tmp_controlled_bit_list = [i for i in controlled_bit_list]
        self.add_GT_ADD_circuit(d+2**self.n_bits-self.N, tmp_controlled_bit_list)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R1)
        tmp_controlled_bit_list.append(self.R1)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R3, self.clean_ancilla, self.dirty_ancilla)
        tmp_controlled_bit_list.remove(self.R1)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R1)
        tmp_controlled_bit_list.append(self.R3)
        self.add_GT_ADD_circuit(self.N, tmp_controlled_bit_list)
        tmp_controlled_bit_list.remove(self.R3)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R1)
        self.add_GT_ADD_circuit(2**self.n_bits-d, tmp_controlled_bit_list)
        tmp_controlled_bit_list.append(self.R1)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R3, self.clean_ancilla, self.dirty_ancilla)
        self.to_clean(self.R3)
        tmp_controlled_bit_list.remove(self.R1)
        self.add_GT_ADD_circuit(d, tmp_controlled_bit_list)
        self.add_ck_not_circuit(tmp_controlled_bit_list, self.R1)
        self.to_clean(self.R1)

    def add_carry_circuit(self, bit_p, ck_plus, ck, tk, controlled_bit_list=[], clean_bit_list=[], dirty_bit_list=[]):
        self.to_dirty(ck_plus)
        if bit_p:
             controlled_bit_list.append(tk)
             self.add_ck_not_circuit(controlled_bit_list, ck_plus, clean_bit_list, dirty_bit_list)
             controlled_bit_list.remove(tk)
             self.add_ck_not_circuit(controlled_bit_list, tk, clean_bit_list, dirty_bit_list)
        
        controlled_bit_list.append(ck)
        controlled_bit_list.append(tk)
        self.add_ck_not_circuit(controlled_bit_list, ck_plus, clean_bit_list, dirty_bit_list)
        controlled_bit_list.remove(ck)
        controlled_bit_list.remove(tk)

    def add_carry_inverse_circuit(self, bit_p, ck_plus, ck, tk, controlled_bit_list=[], clean_bit_list=[], dirty_bit_list=[]):
        controlled_bit_list.append(ck)
        controlled_bit_list.append(tk)
        self.add_ck_not_circuit(controlled_bit_list, ck_plus, clean_bit_list, dirty_bit_list)
        controlled_bit_list.remove(ck)
        controlled_bit_list.remove(tk)

        if bit_p:
            self.add_ck_not_circuit(controlled_bit_list, tk, clean_bit_list, dirty_bit_list)
            controlled_bit_list.append(tk)
            self.add_ck_not_circuit(controlled_bit_list, ck_plus, clean_bit_list, dirty_bit_list)
            controlled_bit_list.remove(tk)
             
        self.to_clean(ck_plus)

    def add_sum_circuit(self, bit_p, ck, tk, controlled_bit_list=[], clean_bit_list=[], dirty_bit_list=[]):
        if bit_p:
            self.add_ck_not_circuit(controlled_bit_list, tk, clean_bit_list, dirty_bit_list)
        controlled_bit_list.append(ck)
        self.add_ck_not_circuit(controlled_bit_list, tk, clean_bit_list, dirty_bit_list)
        controlled_bit_list.remove(ck)

    def add_GT_ADD_circuit(self, p, controlled_bit_list=[]):
        tmp_controlled_bit_list = [i for i in controlled_bit_list]
        for i in range(self.n_bits):
            if p & (1 << i):
                for j in range(i, self.n_bits):
                    tmp_controlled_bit_list.append(self.R2[j])
                self.add_ck_not_circuit(tmp_controlled_bit_list, self.R1, self.clean_ancilla, self.dirty_ancilla)
                for j in range(self.n_bits-1, i-1, -1):
                    tmp_controlled_bit_list.remove(self.R2[j])
                    self.add_ck_not_circuit(tmp_controlled_bit_list, self.R2[j], self.clean_ancilla, self.dirty_ancilla)

    def add_ck_not_circuit(self, controlled_bit_list, target_bit, clean_bit_list=[], dirty_bit_list=[]):
        k = len(controlled_bit_list)
        if k == 0:
            self.add_X_gate(target_bit)
        elif k == 1:
            if controlled_bit_list[0] == target_bit:
                print("error in add_CNOT_gate: using same qubit")
            self.add_c_not_gate(controlled_bit_list[0], target_bit)
        elif k == 2:
            self.add_c2_not_gate(controlled_bit_list[0], controlled_bit_list[1], target_bit)
        else:
            tmp_controlled_bit_list = [i for i in controlled_bit_list]
            if self.use_clean_management:
                tmp_clean_bit_list = [i for i in clean_bit_list]
                tmp_dirty_bit_list = [i for i in dirty_bit_list]
            else:   # cleanビットを使用しない場合
                tmp_clean_bit_list = []
                tmp_dirty_bit_list = [i for i in dirty_bit_list]
                for i in clean_bit_list:
                    tmp_dirty_bit_list.append(i)
            
            if target_bit in tmp_clean_bit_list:
                tmp_clean_bit_list.remove(target_bit)
            if target_bit in tmp_dirty_bit_list:
                tmp_dirty_bit_list.remove(target_bit)
            for i in controlled_bit_list:
                if i in tmp_clean_bit_list:
                    tmp_clean_bit_list.remove(i)
                if i in tmp_dirty_bit_list:
                    tmp_dirty_bit_list.remove(i)
            num_clean_bit = len(tmp_clean_bit_list)

            if num_clean_bit >= 1:
                controlled_bit1 = tmp_controlled_bit_list[0]
                controlled_bit2 = tmp_controlled_bit_list[1]
                using_clean_ancilla = tmp_clean_bit_list[0]
                tmp_controlled_bit_list.remove(controlled_bit1)
                tmp_controlled_bit_list.remove(controlled_bit2)
                self.add_c2_not_gate(controlled_bit1, controlled_bit2, using_clean_ancilla)
                tmp_controlled_bit_list.append(using_clean_ancilla)
                self.add_ck_not_circuit(tmp_controlled_bit_list, target_bit, clean_bit_list, dirty_bit_list)
                self.add_c2_not_gate(controlled_bit1, controlled_bit2, using_clean_ancilla, True)

            else:
                if k - 2 >= len(tmp_dirty_bit_list):
                    print("error in add_ck_not_circuit: lack of dirty ancilla")
                self.add_c2_not_gate(tmp_controlled_bit_list[0], tmp_dirty_bit_list[0], target_bit)
                for i in range(1, k-2):
                    self.add_c2_not_gate(tmp_controlled_bit_list[i], tmp_dirty_bit_list[i], tmp_dirty_bit_list[i-1])
                self.add_c2_not_gate(tmp_controlled_bit_list[k-2], controlled_bit_list[k-1], tmp_dirty_bit_list[k-3])
                for i in range(k-3, 0, -1):
                    self.add_c2_not_gate(tmp_controlled_bit_list[i], tmp_dirty_bit_list[i], tmp_dirty_bit_list[i-1])
                self.add_c2_not_gate(tmp_controlled_bit_list[0], tmp_dirty_bit_list[0], target_bit)

                for i in range(1, k-2):
                    self.add_c2_not_gate(tmp_controlled_bit_list[i], tmp_dirty_bit_list[i], tmp_dirty_bit_list[i-1])
                self.add_c2_not_gate(tmp_controlled_bit_list[k-2], tmp_controlled_bit_list[k-1], tmp_dirty_bit_list[k-3])
                for i in range(k-3, 0, -1):
                    self.add_c2_not_gate(tmp_controlled_bit_list[i], tmp_dirty_bit_list[i], tmp_dirty_bit_list[i-1])
        # self.to_dirty(target_bit)

    def add_c2_not_gate(self, controlled_bit1, controlled_bit2, target_bit, ancilla_status=False):

        if self.decompose_toffoli:
            if controlled_bit1 == target_bit or controlled_bit2 == target_bit or controlled_bit1 == controlled_bit2:
                print("error in add_c2_not_gate: using same qubit")
            self.add_H_gate(target_bit)
            self.add_ck_not_circuit([controlled_bit2], target_bit)
            self.add_Tdag_gate(target_bit)
            self.add_ck_not_circuit([controlled_bit1], target_bit)
            self.add_T_gate(target_bit)
            self.add_ck_not_circuit([controlled_bit2], target_bit)
            self.add_Tdag_gate(target_bit)
            self.add_ck_not_circuit([controlled_bit1], target_bit)
            self.add_T_gate(controlled_bit2)
            self.add_T_gate(target_bit)
            self.add_H_gate(target_bit)
            self.add_ck_not_circuit([controlled_bit1], controlled_bit2)
            self.add_T_gate(controlled_bit1)
            self.add_Tdag_gate(controlled_bit2)
            self.add_ck_not_circuit([controlled_bit1], controlled_bit2)
        else:
            # gate = DenseMatrix(target_bit, [[0, 1],[1, 0]])
            # gate.add_control_qubit(controlled_bit1, 1)
            # gate.add_control_qubit(controlled_bit2, 1)
            #self.circuit.add_gate(gate)
            self.qisc.ccx(controlled_bit1, controlled_bit2, target_bit)
            self.num_c2_not_gate += 1
        if ancilla_status:
            self.to_clean(target_bit)
        else:
            self.to_dirty(target_bit)

    def add_X_gate(self, target_bit, ancilla_status=False):
        #self.circuit.add_X_gate(target_bit)
        self.qisc.x(target_bit)
        if ancilla_status:
            self.to_clean(target_bit)
        else:
            self.to_dirty(target_bit)
        self.num_not_gate += 1

    def add_c_not_gate(self, controlled_bit, target_bit, ancilla_status=False):
        #self.circuit.add_CNOT_gate(controlled_bit, target_bit)
        self.qisc.cx(controlled_bit, target_bit)
        if ancilla_status:
            self.to_clean(target_bit)
        else:
            self.to_dirty(target_bit)
        self.num_c_not_gate += 1
    
    def add_T_gate(self, target_bit, ancilla_status=False):
        #self.circuit.add_T_gate(target_bit)
        self.qisc.t(target_bit)
        if ancilla_status:
            self.to_clean(target_bit)
        else:
            self.to_dirty(target_bit)
        self.num_t_gate += 1

    def add_Tdag_gate(self, target_bit, ancilla_status=False):
        #self.circuit.add_Tdag_gate(target_bit)
        self.qisc.tdg(target_bit)
        if ancilla_status:
            self.to_clean(target_bit)
        else:
            self.to_dirty(target_bit)
        self.num_tdag_gate += 1

    ##############################################################
    #   QFT
    ##############################################################

    def add_qft_circuit(self, n, position=0, swap=True):
        self.qft_rotations(n, position)
        if swap:
            self.swap_registers(n, position)

    def qft_rotations(self, n, position):
        if n == 0:
            return
        n = n - 1

        self.add_H_gate(n+position)
        for i in range(n):
            self.add_CkR_gate(math.pi/2**(n-i), [i+position], n+position) 
        self.qft_rotations(n, position)


    def add_CkR_gate(self, param, controlled_bit_list, target_bit):
        # gate = DenseMatrix(target_bit, [[1,0],[0,exp((0 + 1j) * param)]])
        # for control_bit in controlled_bit_list:
        #     gate.add_control_qubit(control_bit, 1) # 制御ビットが1のときに適用 
        #self.circuit.add_gate(gate)
        self.qisc.mcp(param, controlled_bit_list, target_bit)

        number_of_controlled_bit = len(controlled_bit_list)
        if number_of_controlled_bit == 0:
            self.num_R_gate += 1
        elif number_of_controlled_bit == 1:
            self.num_CR_gate += 1
        elif number_of_controlled_bit == 2:
            self.num_C2R_gate += 1
        else:
            print("Error in add_CkR_gate: number_of_controlled_bit over 2")    

    def swap_registers(self, n, position):
        for i in range(n//2):
            self.add_ck_not_circuit([i+position], n-i-1+position)
            self.add_ck_not_circuit([n-i-1+position], i+position)
            self.add_ck_not_circuit([i+position], n-i-1+position)

    def add_inverse_qft_circuit(self, n, position=0, swap=True):
        if swap:
            self.swap_registers(n, position)
        self.inverse_qft_rotations(-1, n, position)

    def inverse_qft_rotations(self, i, n, position):

        if i == n-1:
            return
        i = i + 1

        for j in range(i-1, -1, -1):
            self.add_CkR_gate(-math.pi/2**(i-j), [j+position], i+position) 

        self.add_H_gate(i+position)
        self.inverse_qft_rotations(i, n, position)

def R_y_matrix(theta):
    return numpy.array([[math.cos(theta/2), math.sin(theta/2)],[-math.sin(theta/2), math.cos(theta/2)]])

def R_z_matrix(alpha):
    return numpy.array([[exp((0 + 1j) * alpha/2), 0], [0, exp((0 + 1j) * (-alpha/2))]])

def Phi_matrix(delta):
    return numpy.array([[exp((0 + 1j) * delta), 0], [0, exp((0 + 1j) * delta)]])
