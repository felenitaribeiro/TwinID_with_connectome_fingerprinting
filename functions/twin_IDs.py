import xlrd

#Index ID python
workbook = xlrd.open_workbook('./../data/List_twin_pairs.xlsx')
sheet = workbook.sheet_by_index(0)

twin_pairs_ID=[]
for rowx in range(sheet.nrows):
    cols = sheet.row_values(rowx)
    twin_pairs_ID.append(cols)


