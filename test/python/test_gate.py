import pytest
from qdd import pyQDD
from qdd.qdd_backend import (
    _qiskit_gates_1q_0param,
    _qiskit_gates_1q_1param,
    _qiskit_gates_1q_2param,
    _qiskit_gates_1q_3param,
    _qiskit_gates_1q_4param,
    _qiskit_gates_2q_0param,
    _qiskit_gates_2q_1param,
)
from qdd.qdd_sampler import Sampler as qdd_sampler
from qiskit import QuantumCircuit
import qiskit.circuit.library as qiskit_gates
from qiskit.primitives import Sampler as qiskit_sampler
import numpy as np
import math
import random
import scipy.stats
import itertools

sampler_qiskit = qiskit_sampler()
sampler_qdd = qdd_sampler(run_options={"shots": None})

qc_size = 6


def test_1q_0param_gates():
    for qis, _ in _qiskit_gates_1q_0param.items():
        if qis == qiskit_gates.MCXGate:
            continue
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        for _ in range(3):
            targets = random.sample(range(qc_size), num_qubits)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qis_gate, targets)
            qc.measure_all()
            job_qis = sampler_qiskit.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_qdd = sampler_qdd.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_result_qis = job_qis.result()
            job_result_qdd = job_qdd.result()
            for dist_qis, dist_qdd in zip(
                job_result_qis.quasi_dists, job_result_qdd.quasi_dists
            ):
                dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_1q_1param_gates():
    for qis, _ in _qiskit_gates_1q_1param.items():
        if qis == qiskit_gates.MCPhaseGate:
            continue
        for _ in range(10):
            para = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = random.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                job_qis = sampler_qiskit.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_qdd = sampler_qdd.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_result_qis = job_qis.result()
                job_result_qdd = job_qdd.result()
                for dist_qis, dist_qdd in zip(
                    job_result_qis.quasi_dists, job_result_qdd.quasi_dists
                ):
                    dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                    dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                    assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_1q_2param_gates():
    for qis, _ in _qiskit_gates_1q_2param.items():
        for _ in range(20):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = random.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                job_qis = sampler_qiskit.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_qdd = sampler_qdd.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_result_qis = job_qis.result()
                job_result_qdd = job_qdd.result()
                for dist_qis, dist_qdd in zip(
                    job_result_qis.quasi_dists, job_result_qdd.quasi_dists
                ):
                    dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                    dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                    assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_1q_3param_gates():
    for qis, _ in _qiskit_gates_1q_3param.items():
        for _ in range(30):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            para3 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = random.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                job_qis = sampler_qiskit.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_qdd = sampler_qdd.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_result_qis = job_qis.result()
                job_result_qdd = job_qdd.result()
                for dist_qis, dist_qdd in zip(
                    job_result_qis.quasi_dists, job_result_qdd.quasi_dists
                ):
                    dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                    dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                    assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_1q_4param_gates():
    for qis, _ in _qiskit_gates_1q_4param.items():
        for _ in range(40):
            para1 = random.uniform(-math.pi, math.pi)
            para2 = random.uniform(-math.pi, math.pi)
            para3 = random.uniform(-math.pi, math.pi)
            para4 = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para1, para2, para3, para4)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = random.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                job_qis = sampler_qiskit.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_qdd = sampler_qdd.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_result_qis = job_qis.result()
                job_result_qdd = job_qdd.result()
                for dist_qis, dist_qdd in zip(
                    job_result_qis.quasi_dists, job_result_qdd.quasi_dists
                ):
                    dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                    dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                    assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_2q_0param_gates():
    for qis, _ in _qiskit_gates_2q_0param.items():
        qis_gate = qis()
        num_qubits = qis_gate.num_qubits
        for _ in range(3):
            targets = random.sample(range(qc_size), num_qubits)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qis_gate, targets)
            qc.measure_all()
            job_qis = sampler_qiskit.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_qdd = sampler_qdd.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_result_qis = job_qis.result()
            job_result_qdd = job_qdd.result()
            for dist_qis, dist_qdd in zip(
                job_result_qis.quasi_dists, job_result_qdd.quasi_dists
            ):
                dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_2q_1param_gates():
    for qis, _ in _qiskit_gates_2q_1param.items():
        for _ in range(10):
            para = random.uniform(-math.pi, math.pi)
            qis_gate = qis(para)
            num_qubits = qis_gate.num_qubits
            for _ in range(3):
                targets = random.sample(range(qc_size), num_qubits)
                qc = QuantumCircuit(qc_size)
                qc.h(range(qc_size))
                qc.append(qis_gate, targets)
                qc.measure_all()
                job_qis = sampler_qiskit.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_qdd = sampler_qdd.run(
                    circuits=[qc], parameter_values=[[]], parameters=[[]]
                )
                job_result_qis = job_qis.result()
                job_result_qdd = job_qdd.result()
                for dist_qis, dist_qdd in zip(
                    job_result_qis.quasi_dists, job_result_qdd.quasi_dists
                ):
                    dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                    dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                    assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)


def test_unitary():
    bit = 3
    for _ in range(10):
        random_matrix = scipy.stats.unitary_group.rvs(2**bit)
        for _ in range(3):
            targets = random.sample(range(qc_size), bit)
            qc = QuantumCircuit(qc_size)
            qc.h(range(qc_size))
            qc.append(qiskit_gates.UnitaryGate(random_matrix), targets)
            qc.measure_all()
            job_qis = sampler_qiskit.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_qdd = sampler_qdd.run(
                circuits=[qc], parameter_values=[[]], parameters=[[]]
            )
            job_result_qis = job_qis.result()
            job_result_qdd = job_qdd.result()
            for dist_qis, dist_qdd in zip(
                job_result_qis.quasi_dists, job_result_qdd.quasi_dists
            ):
                dist_qis = {k: v for k, v in dist_qis.items() if v != 0}
                dist_qdd = {k: v for k, v in dist_qdd.items() if v != 0}
                assert dist_qis == pytest.approx(dist_qdd, rel=1e-6)
