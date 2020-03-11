def listToString(s):  
    
    # initialize an empty string 
    str1 = ""  
    
    # traverse in the string   
    for ele in s:  
        str1 += ele   
    
    # return string   
    return str1
    
def solution(n, b):
    #Your code here
    k = len(n)
    jobidasc = sorted(n)
    jobiddesc = sorted(n, reverse = True)
    
    y = listToString(jobidasc)
    x = listToString(jobiddesc)
    print('Ascending ' + x)
    print('Descending ' + y)

    z = int(x)-int(y)
    z = str(z)
    z = z.zfill(k)
    print(z)

solution('12889', 12)