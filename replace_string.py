#!/usr/bin/env python3

"""

"""

import sys, math, string, re 

#########################################################################################################

tag_list = []


tag_list.append( (  'subroutine','SUBROUTINE' )  )
tag_list.append( (  '^ *function','FUNCTION' )  )
tag_list.append( (  '_FUNCTION','function' )  )
tag_list.append( (  'return','RETURN' )  )
tag_list.append( (  'enddo','END DO' )  )
tag_list.append( (  'end do','END DO' )  )
tag_list.append( (  'call','CALL' )  )
tag_list.append( (  'iCALL','icall' )  )
#tag_list.append( (  '^stop','STOP' )  )
tag_list.append( (  'stop','STOP' )  )
tag_list.append( (  ' stop',' STOP' )  )
tag_list.append( (  'do *$','DO ' )  )
tag_list.append( (  ' do ',' DO ' )  )
tag_list.append( (  r"dabs\(", r"DABS (" )  )
tag_list.append( (  r"any\(", r"ANY (" )  )
tag_list.append( (  r"abs\(", r"ABS (" )  )
tag_list.append( (  r"abs \(", r"ABS (" )  )
tag_list.append( (  r"real\(",r"REAL (" )  )
tag_list.append( (  r"real \(",r"REAL (" )  )
tag_list.append( (  r"float\(",r"FLOAT (" )  )
tag_list.append( (  r"integer\(",r"INTEGER (" )  )
tag_list.append( (  r"nint\(",r"NINT (" )  )
tag_list.append( (  r"int\(",r"INT (" )  )
tag_list.append( (  'integer','INTEGER' )  )
tag_list.append( (  r"kind\(",r"KIND (" )  )
tag_list.append( (  'kind=','KIND=' )  )
#tag_list.append( (  '^use ','USE ' )  )
tag_list.append( (  'use ','USE ' )  )
tag_list.append( (  ' use ',' USE ' )  )
tag_list.append( (  r"\)then",r") THEN" )  )
tag_list.append( (  'then','THEN' )  )
tag_list.append( (  'implicit','IMPLICIT' )  )
tag_list.append( (  'intent','INTENT' )  )
tag_list.append( (  r"\(in\)",r"(IN)" )  )
tag_list.append( (  r"\(out\)", r"(OUT)" )  )
tag_list.append( (  'dimension','DIMENSION' )  )
tag_list.append( (  r"elseif\(", r"ELSE IF (" )  )
tag_list.append( (  'else','ELSE' )  )
tag_list.append( (  'end if','END IF' )  )
tag_list.append( (  'endif','END IF' )  )
#tag_list.append( (  'ENDif','END IF' )  )
#tag_list.append( (  'END if','END IF' )  )
#tag_list.append( (  '^if(','IF (' )  )
#tag_list.append( (  r"^ *if\(", r"IF \(" )  )
tag_list.append( (  r"if\(", r"IF (" )  )
tag_list.append( (  r" if\(", r" IF (" )  )
tag_list.append( (  r" if \(", r" IF (" )  )
tag_list.append( (  r"sqrt\(", r"SQRT (" )  )
tag_list.append( (  r"sqrt \(", r"SQRT (" )  )
tag_list.append( (  r"^interface",r"INTERFACE" )  )
tag_list.append( (  r" interface",r" INTERFACE" )  )
tag_list.append( (  'procedure','PROCEDURE' )  )
tag_list.append( (  'pointer','POINTER' )  )
tag_list.append( (  'nopass','NOPASS' )  )
tag_list.append( (  ',parameter',',PARAMETER' )  )
tag_list.append( (  '_PARAMETER','_parameter' )  )
tag_list.append( (  ',parameter,',',PARAMETER,' )  )
tag_list.append( (  r"^type", r"TYPE" )  )
tag_list.append( (  r" type", r" TYPE" )  )
tag_list.append( (  'contains','CONTAINS' )  )
tag_list.append( (  r"isnan\(", r"ISNAN (" )  )
tag_list.append( (  r"exp\(", r"EXP (" )  )
tag_list.append( (  r"min\(", r"MIN (" )  )
tag_list.append( (  r"max\(", r"MAX (" )  )
tag_list.append( (  r"min \(", r"MIN (" )  )
tag_list.append( (  r"max \(", r"MAX (" )  )
#tag_list.append( (  '^logical','LOGICAL' )  )
tag_list.append( (  'logical','LOGICAL' )  )
tag_list.append( (  ' logical',' LOGICAL' )  )
tag_list.append( (  r"mod\(", r"MOD (" )  )
#tag_list.append( (  '^data','DATA' )  )
tag_list.append( (  'data','DATA' )  )
tag_list.append( (  ' data',' DATA' )  )
tag_list.append( (  'DATA_','data_' )  )
tag_list.append( (  r"acos\(", r"ACOS (" )  )
tag_list.append( (  r"cos\(", r"COS (" )  )
tag_list.append( (  r"sin\(", r"SIN (" )  )
tag_list.append( (  'recursive','RECURSIVE' )  )
tag_list.append( (  r"result \(", r"RESULT (" )  )
tag_list.append( (  r"result\(", r"RESULT (" )  )
tag_list.append( (  r"associated\(", r"ASSOCIATED (" )  )
tag_list.append( (  r"deallocate\(", r"DEALLOCATE (" )  )
tag_list.append( (  r"allocate\(", r"ALLOCATE (" )  )
tag_list.append( (  r"allocated\(", r"ALLOCATED (" )  )
tag_list.append( (  r"allocatable", r"ALLOCATABLE" )  )
#tag_list.append( (  "allocate\\(", "ALLOCATE \\(" )  )
#tag_list.append( (  "allocated\\(", "ALLOCATED \\(" )  )
tag_list.append( (  r"select case\(", r"SELECT CASE (" )  )
tag_list.append( (  r"end *select",r"END SELECT" )  )
tag_list.append( (  'case default','CASE DEFAULT' )  )
tag_list.append( (  r"case\(", r"CASE (" )  )
tag_list.append( (  r"case \(", r"CASE (" )  )
tag_list.append( (  r"null\(", r"NULL (" )  )
tag_list.append( (  'public','PUBLIC' )  )
tag_list.append( (  r"extends\(", r"EXTENDS (" )  )
tag_list.append( (  r"class\(", r"CLASS (" )  )
tag_list.append( (  r"\(inout\)", r"(INOUT)" )  )
tag_list.append( (  r"write\(", r"WRITE (" )  )
tag_list.append( (  r"read\(", r"READ (" )  )
tag_list.append( (  r"open\(", r"OPEN (" )  )
tag_list.append( (  r"close\(", r"CLOSE (" )  )
tag_list.append( (  r"rewind\(", r"REWIND (" )  )
tag_list.append( (  r"character\(", r"CHARACTER (" )  )
tag_list.append( (  'len=','LEN=' )  )
tag_list.append( (  'iostat','IOSTAT' )  )
tag_list.append( (  'cycle','CYCLE' )  )
tag_list.append( (  'continue','CONTINUE' )  )
tag_list.append( (  r"index\(", r"INDEX (" )  )
tag_list.append( (  r"trim\(", r"TRIM (" )  )
tag_list.append( (  r"system\(", r"SYSTEM (" )  )
#tag_list.append( (  '^ *module','MODULE' )  )
tag_list.append( (  r"^ *module", r"MODULE" )  )
tag_list.append( (  r" module", r" MODULE" )  )
tag_list.append( (  'mpi_module','mpi_module' )  )
tag_list.append( (  'abstract','ABSTRACT' )  )
tag_list.append( (  'target','TARGET' )  )
tag_list.append( (  'select type','SELECT TYPE' )  )
tag_list.append( (  'type is','TYPE IS' )  )
tag_list.append( (  'random_number','RANDOM_NUMBER' )  )
tag_list.append( (  'Random_Number','RANDOM_NUMBER' )  )
tag_list.append( (  'do while','DO WHILE' )  )
tag_list.append( (  r"^external",r"EXTERNAL" )  )
tag_list.append( (  r" external",r" EXTERNAL" )  )
tag_list.append( (  ' end',' END' )  )
tag_list.append( (  '^end *','END ' )  )
tag_list.append( (  r"lmdif\(", r"lmdif(" )  )



#########################################################################################################


inputfilename = sys.argv[1]
inputfile = open( inputfilename , 'r' )


contents = inputfile.readlines()
inputfile.close()



outputfilename =  'UC_' + inputfilename 

outputfile = open( outputfilename, 'w' )

#########################################################################################################




j = 0 

while  j < len(contents) :


    line = contents[j]

    #print( 'j, line ', j, line[:-1] ) 

    if  re.match( '^ *!', line ) : 
        #print( 'j, COMMENT ', j  ) 
        outputfile.write( line ) 
        j = j + 1   
        continue

    
        
    for item in tag_list :

        #print( 'item ', item ) 

        #if item[0] in line:
        #if re.search( item[0] ,  line) :
        pat = item[0]
        if re.search( pat ,  line, re.IGNORECASE ) :


            #print( "item[0], item[1] ", item[0], item[1] ) 

            line = re.sub( pat, item[1], line, re.IGNORECASE  ) 

            #re.sub( item[0], item[1], line ) 

            #line = line.replace(  item[0], item[1] ) 

            #break
                

    
    #print( 'j, line ', j, line[:-1] ) 
    outputfile.write( line ) 

    j = j + 1   



outputfile.close()



