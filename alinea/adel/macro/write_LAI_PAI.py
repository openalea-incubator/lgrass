def Write_LAI_PAI(Tsum,Gc,CAI,CAIg,LAI,LAIg,filename):
    '''
    returns a string representation of the content of the table. 
    '''

    fout = open(filename,"a")
    
    fout.write(str(Tsum))
    fout.write(' , ')
    fout.write(str(Gc))
    fout.write(' , ')
    fout.write(str(CAI))
    fout.write(' , ')
    fout.write(str(CAIg))
    fout.write(' , ')
    fout.write(str(LAI))
    fout.write(' , ')
    fout.write(str(LAIg))
    fout.write('\n')
    fout.close()
    
    return filename,
