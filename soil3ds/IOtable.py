from os.path import join

def N_lignes (fichier) :
    """compte le nombre de lignes d'un fichier (compte le nombre d'elements de la liste readlines()"""

    position_ini = fichier.tell()
    fichier.seek(0)
    N_tot_lignes = len (fichier.readlines())
    fichier.seek(position_ini)
    return N_tot_lignes

def transcript_csv_str (fichier) :
     """ transcrit une ligne d'un fichier.csv en liste de float """
     ligne_ch = fichier.readline()
     liste, i, caract = [], 0, '' #initialise une liste vide et un indice a 0 et une chaine de caractere vide 

     if ligne_ch == '' : #verifie si la chaine est vide
         return 'chaine vide'
   
     while ligne_ch[i] !=  '\n': #boucle pour une valeur en .csv (separateur ;)
         ch = ligne_ch[i]
         if ligne_ch[i] == ',' : #converti separateur decimal en point
            ch = '.'
         if ligne_ch[i] == ';'  :
            liste.append (caract)
            caract = ''
         if ligne_ch[i] != ';' :
            caract = caract + ch
         i=i+1

     liste.append (caract) #ajoute la derniere valeur
     return liste

def table_csv_str (fichier) :
    """transcrit un fichier .csv de nombres en un tableau (liste de liste) de nombre reels """

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)
    
    for i in range(n) :
        liste.append (transcript_csv_str (fichier))

    return liste

def transcript_csv (fichier) :
     """ transcrit une ligne d'un fichier.csv en liste de float """
     ligne_ch = fichier.readline()
     liste, i, caract = [], 0, '' #initialise une liste vide et un indice a 0 et une chaine de caractere vide

     if ligne_ch == '' : #verifie si la chaine est vide
         return 'chaine vide'
   
     while ligne_ch[i] !=  '\n': #boucle pour une valeur en .csv (separateur ;)
         ch = ligne_ch[i]
         if ligne_ch[i] == ',' : #converti separateur decimal en point
            ch = '.'
         if ligne_ch[i] == ';'  :
            liste.append (float(caract))
            caract = ''
         if ligne_ch[i] != ';' :
            caract = caract + ch
         i=i+1

     liste.append (float(caract)) #ajoute la derniere valeur

     return liste

def table_csv (fichier) :
    """transcrit un fichier .csv de nombres en un tableau (liste de liste) de nombre reels """

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)
    
    for i in range(n) :
        liste.append (transcript_csv (fichier))

    return liste

def transcript_txt (fichier) :
     """ transcrit une ligne d'un fichier.txt en liste de chaine """
     ligne_ch = fichier.readline()
     liste, i, caract = [], 0, '' #initialise une liste vide et un indice a 0 et une chaine de caractere vide 

     if ligne_ch == '' : #verifie si la chaine est vide
         return 'chaine vide'
   
     while ligne_ch[i] !=  '\n': #boucle pour une valeur en .csv (separateur ;)
         ch = ligne_ch[i]
         if (ligne_ch[i] == ' ' or ligne_ch[i] == "\t") and caract != '' :
            liste.append (caract)
            caract = ''
         if ligne_ch[i] != ' ' and ligne_ch[i] != "\t":
            caract = caract + ch
         i=i+1

     liste.append (caract) #ajoute la derniere valeur
     return liste


def table_txt (fichier) :
    """transcrit un fichier .txt un tableau (liste de liste) de chaines """

    fichier.seek(0)
    liste, n = [], N_lignes(fichier)
    
    for i in range(n) :
        liste.append (transcript_txt (fichier))

    return liste

def ecriture_csv (table, fichier) :
    """ ecrit table (liste de liste) de chiffre dans un fichier csv """
    for i in range (len(table)):
        for j in range (len(table[i][0:])-1) :
            fichier.write(str(table[i][j]))
            fichier.write(';')
        fichier.write(str(table[i][j+1]))
        fichier.write('\n')
    fichier.close()

	
def ecriture_csv_fromlist (ltable, fichier) :
    """ ecrit LISTE de table (LISTE de liste de liste) de chiffre dans un fichier csv """
    for n in range (len(ltable)):
        for i in range (int(n!=0),len(ltable[n])):#int(n!=0)permet deviter de recopier la ligne de header a chaque fois
            fichier.write(str(n))
            fichier.write(';')
            for j in range (len(ltable[n][i][0:])-1) :
                fichier.write(str(ltable[n][i][j]))
                fichier.write(';')
            fichier.write(str(ltable[n][i][j+1]))
            fichier.write('\n')
    fichier.close()



def ecriture_txt (table, fichier) :
    """ ecrit table (liste de liste) de chiffre dans un fichier txt """
    for i in range (len(table)):
        for j in range (len(table[i][0:])-1) :
            fichier.write(str(table[i][j]))
            fichier.write(' ')
        fichier.write(str(table[i][j+1]))
        fichier.write('\n')
    fichier.close()

def copie_partielle (fichier,out,n,m):
    """ recopie de fichier des lignes n a m dans le fichier out """
    for i in range(m):
        ligne_ch = fichier.readline()
        if i>=n:
            out.writelines(ligne_ch)


def conv_dataframe(tab):
    """ converti liste de liste en dictionnaire; prend la cle comme le pemier element de la liste"""
    """ format a priori compatible pour conversion en data.frame R"""
    dat = {}
    for i in range(len(tab)):
        dat[str(tab[i][0])] = tab[i][1:]

    return dat #r.as_data_frame(dat)

def conv_list(tab):
    """ converti dictionnaireen liste de liste en ;  cle comme pemier element de la liste"""
    """ format compatible pour mes_csv"""
    dat = []
    for i in list(tab.keys()):
        v = [i]
        dat.append(v)

    count = 0
    for i in list(tab.keys()):
        for j in range(len(tab[i])):
            dat[count].append(tab[i][j])

        count = count+1

    return dat 

def extract_dataframe(dat, ls_cles, cle, val=None, oper='egal'):
    """ extrait dans listes de cles ls_cles les lignes pour lesquelles cle=val; toutes si val=None  
    option pour oper: egal / sup / inf/ supeg / infeg / diff"""
    #cree liste d'index ou cle = val

    id = []
    for i in range(len(dat[cle])):
        if val == None:
            id.append(i)
        else:
            if oper =='egal':
                if dat[cle][i] == val:
                    id.append(i)
            elif oper =='inf':
                if dat[cle][i] < val:
                    id.append(i)
            elif oper =='sup':
                if dat[cle][i] > val:
                    id.append(i)
            elif oper =='infeg':
                if dat[cle][i] <= val:
                    id.append(i)
            elif oper =='supeg':
                if dat[cle][i] >= val:
                    id.append(i)
            elif oper =='diff':
                if dat[cle][i] != val:
                    id.append(i)

    x = {}
    for k in ls_cles: # recupere les paires interessantes
        v = []
        for i in id: # les id respectant cle=val
            v.append(dat[k][i])

        x[k] = v

    return x
    #extract_dataframe(dat, cles, 'geno', geno)
    #extract_dataframe(dat, cles, 'geno')

def t_list(tab):
    """transpose tab"""
    res = []
    for j in range(len(tab[0])):
        v = []
        for i in range(len(tab)):
            v.append(tab[i][j])
        
        res.append(v)

    return res


def conv_list2(tab):
    """ converti dictionnaire avec 1 seule valeur par cle en liste de liste ;  cle comme pemier element de la liste"""
    """ format compatible pour mes_csv"""
    dat = []
    for i in list(tab.keys()):
        v = [i, tab[i]]
        dat.append(v)
    
    return dat 


def write_dict(dict, directory, name):

    try:
        tab = t_list(conv_list(dict))
    except:
        tab = conv_list2(dict)

    out = open(join(directory, name), 'w')#file(join(directory, name), 'w')
    ecriture_csv (tab, out)  
    out.close()
    return join(directory, name)


def write_dicttables(path_file, dic, keys_):
    """ ecrit dans un meme fichiers, differentes tables (liste de liste) valeurs d'un dico pour la liste de cle keys_ """
    f = open(path_file, 'w')#file (path_file, 'w')
    ecriture_csv(dic[keys_[0]], f)
    f.close()
    if len(keys_)>1:
        for i in range (1,len(keys_)):
            f = open(path_file, 'a')#file (path_file, 'a')
            ecriture_csv(dic[keys_[i]], f)
            f.close()

