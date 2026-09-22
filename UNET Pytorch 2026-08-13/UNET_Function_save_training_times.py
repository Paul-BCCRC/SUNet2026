from openpyxl import Workbook

def save_training_times(filename, totaltime, readtime, traintime, predicttime, writetime):
    
    workbook = Workbook()
    sheet = workbook.active
    
    sheet["A1"] = "Total"
    sheet["A2"] = "Read"
    sheet["A3"] = "Train"
    sheet["A4"] = "Predict"
    sheet["A5"] = "Write"
    
    sheet["B1"] = totaltime
    sheet["B2"] = readtime
    sheet["B3"] = traintime
    sheet["B4"] = predicttime
    sheet["B5"] = writetime

    workbook.save(filename)
    
    
    return []