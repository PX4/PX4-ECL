#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 14 12:47:24 2020

@author: roman
"""
from sympy import ccode
from sympy.codegen.ast import float32, real

class CodeGenerator:
    def __init__(self, file_name):
        self.file_name = file_name
        self.file = open(self.file_name, 'w', buffering=0)
        
    def get_ccode(self, expression):
        return ccode(expression, type_aliases={real:float32})
        
    def write_subexpressions(self,subexpressions):
        write_string = ""
        for item in subexpressions:
            write_string = write_string + "float " + str(item[0]) + " = " + self.get_ccode(item[1]) + ";\n"
        
        write_string = write_string + "\n\n"
        self.file.write(write_string)
        
    def write_matrix(self, matrix, identifier, is_symmetric=False):
        write_string = ""
        
        if matrix.shape[0] == 1 or matrix.shape[1] == 1:
            for i in range(0,len(matrix)):
                write_string = write_string + identifier + "(" + str(i) + ") = " + self.get_ccode(matrix[i]) + ";\n"
        else:
            #write_string = "float " + identifier + "[" + str(matrix.shape[0]) + "][" + str(matrix.shape[1]) + "] = {};\n"
            for j in range(0, matrix.shape[1]):
                for i in range(0, matrix.shape[0]):
                    if j >= i or not is_symmetric:
                        write_string = write_string + identifier + "(" + str(i) + "," + str(j) + ") = " + self.get_ccode(matrix[i,j]) + ";\n"

        write_string = write_string + "\n\n"
        self.file.write(write_string)
        
    def close(self):
        self.file.close()