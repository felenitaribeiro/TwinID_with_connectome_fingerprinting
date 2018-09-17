#!/usr/bin/env python3
"""
Script to transform a list of twin identities on excel into a list on Python
"""
import xlrd

# Index ID python
workbook = xlrd.open_workbook('/home/hardleste/Desktop/HCP_DATA_12_2017/python_scripts_ESTUDO1/List_twin_pairs.xlsx')
sheet = workbook.sheet_by_index(0)

twin_pairs_ID = []
for rowx in range(sheet.nrows):
    cols = sheet.row_values(rowx)
    twin_pairs_ID.append(cols)
