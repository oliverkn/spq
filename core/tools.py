import numpy as np
import qutip as qu

def isomorphism():
    pass

def super_isomorphism(isomorphism):
    def eps_phi(op):
        return isomorphism * op * isomorphism.dag()

    return eps_phi


class Factorization:
    def __init__(self):
        pass

    def get_state_isomorphism(self):
        pass

    def get_operator_isomorphism(self):
        pass

