"""
List of tensor operations.
"""

import sympy as sp
import numpy as np

def covariant_derivative(tensor, spacetime, index_type='lower'):
    """
    Calculate the covariant derivative of a tensor.
    The covariant derivative of a tensor T_{ij} is given by:
       ∇_k T_{ij} = ∂_k T_{ij} - Γ^m_{ki} T_{mj} - Γ^m_{kj} T_{im}
    where Γ^m_{ki} are the Christoffel symbols of the second kind.
    """
    Gamma = spacetime.christoffel_symbols()
    coords = spacetime.coords
    n = len(coords)
    
    if index_type == 'lower':
        # T_{ij}
        cov_derivative = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    cov_derivative[k][i][j] = sp.diff(tensor[i, j], coords[k]) - sum(
                        Gamma[m][k][i] * tensor[m, j] + Gamma[m][k][j] * tensor[i, m] for m in range(n)
                    )
    elif index_type == 'upper':
        # T^{ij}
        cov_derivative = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    cov_derivative[k][i][j] = sp.diff(tensor[i, j], coords[k]) + sum(
                        Gamma[i][k][m] * tensor[m, j] + Gamma[j][k][m] * tensor[i, m] for m in range(n)
                    )
    else:
        raise ValueError("index_type must be 'lower' or 'upper'")
    
    return cov_derivative

def lie_derivative(tensor, vector_field, spacetime, index_type='lower'):
    """
    Calculate the Lie derivative of a tensor along a vector field.
    The Lie derivative of a tensor T_{ij} along a vector field X^k is given by:
       L_X T_{ij} = X^k ∂_k T_{ij} + T_{kj} ∂_i X^k + T_{ik} ∂_j X^k
    where X^k is the vector field and T_{ij} is the tensor.
    """
    coords = spacetime.coords
    n = len(coords)
    
    if index_type == 'lower':
        # T_{ij}
        lie_derivative = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                lie_derivative[i][j] = sum(
                    vector_field[k] * sp.diff(tensor[i, j], coords[k]) + tensor[k, j] * sp.diff(vector_field[i], coords[k]) + tensor[i, k] * sp.diff(vector_field[j], coords[k])
                    for k in range(n)
                )
    elif index_type == 'upper':
        # T^{ij}
        lie_derivative = [[[0 for _ in range(n)] for _ in range(n)] for _ in range(n)]
        for i in range(n):
            for j in range(n):
                lie_derivative[i][j] = sum(
                    vector_field[k] * sp.diff(tensor[i, j], coords[k]) - tensor[k, j] * sp.diff(vector_field[i], coords[k]) - tensor[i, k] * sp.diff(vector_field[j], coords[k])
                    for k in range(n)
                )
    else:
        raise ValueError("index_type must be 'lower' or 'upper'")
    
    return lie_derivative

def tensor_contraction(tensor, index1, index2):
    """
    Contract a tensor over two indices.
    The contraction of a tensor T_{ij} over indices i and j is given by:
       T^i_i = g^{ij} T_{ij}
    where g^{ij} is the inverse metric tensor.
    """
    n = tensor.shape[0]
    contracted_tensor = 0
    for i in range(n):
        contracted_tensor += tensor[i, i]
    return contracted_tensor

def tensor_product(tensor1, tensor2):
    """
    Calculate the tensor product of two tensors.
    The tensor product of two tensors T_{ij} and S_{kl} is given by:
       (T ⊗ S)_{ijkl} = T_{ij} S_{kl}
    where T_{ij} and S_{kl} are the components of the tensors.
    """
    shape1 = tensor1.shape
    shape2 = tensor2.shape
    product_shape = shape1 + shape2
    product_tensor = np.zeros(product_shape, dtype=object)
    
    for i in range(shape1[0]):
        for j in range(shape1[1]):
            for k in range(shape2[0]):
                for l in range(shape2[1]):
                    product_tensor[i, j, k, l] = tensor1[i, j] * tensor2[k, l]
    
    return product_tensor

def symmetrize_tensor(tensor, indices):
    """
    Symmetrize a tensor over specified indices.
    The symmetrization of a tensor T_{ij} over indices i and j is given by:
       T_{(ij)} = (1/2) * (T_{ij} + T_{ji})
    where T_{ij} is the original tensor and T_{(ij)} is the symmetrized tensor.
    """
    n = tensor.shape[0]
    symmetrized_tensor = np.zeros(tensor.shape, dtype=object)
    
    for i in range(n):
        for j in range(n):
            symmetrized_tensor[i, j] = (tensor[i, j] + tensor[j, i]) / 2
    
    return symmetrized_tensor

def antisymmetrize_tensor(tensor, indices):
    """
    Antisymmetrize a tensor over specified indices.
    The antisymmetrization of a tensor T_{ij} over indices i and j is given by:
       T_{[ij]} = (1/2) * (T_{ij} - T_{ji})
    where T_{ij} is the original tensor and T_{[ij]} is the antisymmetrized tensor.
    """
    n = tensor.shape[0]
    antisymmetrized_tensor = np.zeros(tensor.shape, dtype=object)
    
    for i in range(n):
        for j in range(n):
            antisymmetrized_tensor[i, j] = (tensor[i, j] - tensor[j, i]) / 2
    
    return antisymmetrized_tensor

def trace(tensor):
    """
    Calculate the trace of a tensor.
    The trace of a tensor T_{ij} is given by:
       Tr(T) = T^i_i = g^{ij} T_{ij}
    where g^{ij} is the inverse metric tensor and T_{ij} is the tensor.
    """
    n = tensor.shape[0]
    trace_value = 0
    for i in range(n):
        trace_value += tensor[i, i]
    return trace_value

def determinant(tensor):
    """
    Calculate the determinant of a tensor.
    The determinant of a tensor T_{ij} is given by:
       det(T) = ε^{i_1 i_2 ... i_n} T_{i_1 j_1} T_{i_2 j_2} ... T_{i_n j_n}
    where ε^{i_1 i_2 ... i_n} is the Levi-Civita symbol and T_{ij} is the tensor.
    """
    return sp.simplify(tensor.det())

def inverse(tensor):
    """
    Calculate the inverse of a tensor.
    The inverse of a tensor T_{ij} is given by:
       T^{-1}_{ij} = (1/det(T)) * adj(T)_{ij}
    where det(T) is the determinant of the tensor and adj(T)_{ij} is the adjugate of the tensor.
    """
    return sp.simplify(tensor.inv())

def transpose(tensor):
    """
    Calculate the transpose of a tensor.
    The transpose of a tensor T_{ij} is given by:
       T^T_{ij} = T_{ji}
    where T_{ij} is the original tensor and T^T_{ij} is the transposed tensor.
    """
    return sp.simplify(tensor.T)