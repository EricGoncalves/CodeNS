#!/usr/bin/env python2
# -*- coding: utf8 -*-

from os import listdir,path
from os.path import isfile, join
import sys

import fileinput

import subprocess
from subprocess import Popen, check_output,call,CalledProcessError,PIPE

import re
## Import graphviz
#import gv
## Import pygraph
#from pygraph.classes.digraph import digraph
#from pygraph.readwrite.dot import write
#from pygraph.algorithms.cycles import find_cycle
#from pygraph.algorithms.critical import transitive_edges
#from pygraph.algorithms.accessibility import *
#from pygraph.algorithms.sorting import topological_sorting

#implicit none everywhere

def matched(string, regexp):
  if regexp and string:
    return re.search(regexp, string) is not None
  return None

#implicit=False
#module=False

#listfiles = [ f for f in listdir(".") if isfile(join(".",f)) ]
#for fich in listfiles:
#  fileExtension = path.splitext(fich)[1]
#  if fileExtension == ".f90" :
#    if implicit:
#      gfort=Popen(['gfortran', "-O0","-fimplicit-none","-fmax-errors=500","-c",fich],stderr=subprocess.STDOUT,stdout=PIPE)
#      output= gfort.communicate()[0]
#      if gfort.returncode != 0:
#        i=0
#        for line in fileinput.input(fich):
#          if matched(line,"subroutine"):
#            i= fileinput.filelineno()
#          if matched(line," use "):
#            i= fileinput.filelineno()
#          if matched(line,"function"):
#            i= fileinput.filelineno()
#        for line in fileinput.input(fich,inplace=1):
#          print line,
#          if fileinput.filelineno()==i:
#            line2=""
#            line3=""
#            for line1 in output.splitlines():
#              if matched(line1,"IMPLICIT"):
#                var=line1.split("'")[1]
#                if var[0] in "ijklmn":
#                  print "integer :: "+var
#                else:
#                  print "double precision :: "+var
#              if matched(line1,"Variable type is UNKNOWN"):
#                var=line3.split("(")[0].replace(" ","")
#                if var[0] in "ijklmn":
#                  print "integer :: "+var
#                else:
#                  print "double precision :: "+var
#              line3=line2
#              line2=line1
#    if module:
#      modules=""
#      for line in fileinput.input(fich):
#        if matched(line,"^ *subroutine"):
#          i= fileinput.filelineno()
#        if matched(line,"^ *use "):
#          i= fileinput.filelineno()
#        if matched(line,"^ *function"):
#          i= fileinput.filelineno()
#        if not matched(line,"^ *!"):
#          if matched(line,"call "):
#            modules=modules+" "+line.split("(")[0].split("call ")[1].replace(" ","")
#      if modules:
#        for line in fileinput.input(fich,inplace=1):
#          print line,
#          if fileinput.filelineno()==i:
#            for mod in set(sorted(modules.split(" "))):
#              if mod:
#                print "use mod_"+mod

#line2=""
#line3=""
#line4=""
#line5=""
#for line1 in open("compil.log"):
#  if "Unused dummy argument" in line1:
#    var=line1.split("'")[1]
#    fich=line5.split(":")[0]
#    test=False
#    test1=False
#    for line6 in open(fich):
#      if "dimension" in line6:
#        test1=True
#      if test1:
#        if var in line6:
#          test=True
#      if not "&" in line6:
#        test1=False
#    if not test:
#      print fich+" : "+str(var)
#  line5=line4
#  line4=line3
#  line3=line2
#  line2=line1



# dot generator

#def clean(chaine):
#  return chaine.split("\n")[0].strip().replace('mod_','')

#def clean_list(monster):
#    monster=list(set(monster))
#    newmonster=[]
#    for elem in monster:
##      if not "modules" in elem:
#        newmonster.append(elem)
#    return newmonster

#def depgraph(fich):
#    gr = digraph()
#    gr.add_nodes([fich])
#    routines=[]
#    for line in open(fich+".f90"):
#      if not matched(line,"^ *!"):
#        if "use " in line:
#          routines.append(clean(line.split("use ")[1]))
#    routines2=[]
#    for routine in list(set(routines)):
#      if not isfile(routine+".f90"):
#        routine="modules"
#      routines2.append(routine)
#    for routine in list(set(routines2)):
#      if not gr.has_node(routine):
#        newgr=depgraph(routine)
#        gr.add_graph(newgr)
#      if gr.has_node(routine):
#        gr.add_edge((fich, routine))
#    return gr

## Graph creation
#gr1 = depgraph("solve")

##cycles=find_cycle(gr1)
##if cycles:
##  print "cycles : "
##  print cycles
#gr2=gr1
##else:
##  tr_edges=transitive_edges(gr1)
##  tr_edges=clean_list(tr_edges)
##  if tr_edges:
##    print "tr_edges : "
##    print tr_edges
##    gr2 = depgraph(tr_edges[0][0])
##  else:
##    gr2=gr1

#gr2.del_node("modules")
#gr2.del_node("solve")
## Draw as PS
#dot = write(gr2)
##gvv = gv.readstring(dot)
##gv.layout(gvv,'dot')
##gv.render(gvv,'ps','graph.ps')
#print dot

#print "cmake_minimum_required(VERSION 2.6)"
#print "project(CodeNS Fortran)"
#for src in gr2.nodes():
#  if src=="solve":
#    print "add_executable(solver_air "+src+".f90)"
#  else:
#    print "add_library("+src+" SHARED "+src+".f90)"
#for dep in gr2.edges():
#  if dep[0]=="solve":
#    src="solver_air"
#  else:
#    src=dep[0]
#  print "target_link_libraries("+src+" "+dep[1]+")"



##print "target_link_libraries(solve lib3)
##print "file(GLOB source_files *.f90)"
##print "add_executable(solver_air ${source_files})"


#def pprint(prepend,liste):
#  list1=sorted(liste)
#  larg=len(max(list1, key=len))
#  for index, elem in enumerate(list1):
#    if index % 5 ==0:
#      sys.stdout.write(prepend)
#    else:
#      sys.stdout.write(",")
#    sys.stdout.write(elem.rjust(larg))
#    if index % 5 ==4:
#      sys.stdout.write("\n")
#  if not index % 5 ==4:
#    print ""



#def get_in_list(liste,typevar,allocatable,replace,args):
#  result=[]
#  for key, value in liste.iteritems():
#    if value[0]==typevar:
#      var=key
#      allocatablize=False
#      if value[1]:
#        var=var+"("
#        if replace and not key in args:
#          for dim in value[1]:
#            for dim1 in dim.split("*"):
#              if dim1 in excluded_dim:
#                allocatablize=True
#        for index, dim in enumerate(value[1]):
#          if allocatablize:
#            dim=":"
#          if index==0:
#            var=var+dim
#          else:
#            var=var+","+dim
#        var=var+")"
#      if ":" in var and allocatable:
#        result.append(var)
#      if not ":" in var and not allocatable:
#        result.append(var)
#  return result

#def count_units(fich):
#  j=0
#  for line in open(fich):
#    if matched(line,"^ *subroutine") or \
#       matched(line,"^ *module")     or \
#       matched(line,"^ *function")   or \
#       matched(line,"^ *program")  :
#       j=j+1
#  return j

#def check_start(line,i):
#  global k,read,inarg
#  if matched(line,"^ *subroutine") or \
#     matched(line,"^ *module")     or \
#     matched(line,"^ *function")   or \
#     matched(line,"^ *program")  :
#     k=k+1
#     if(k==i):
#       read=True
#       if matched(line,"^ *subroutine") or \
#          matched(line,"^ *function") :
#         inarg=True
#  if matched(line,"^ *end subroutine") or \
#     matched(line,"^ *end module")     or \
#     matched(line,"^ *contains")       or \
#     matched(line,"^ *end function")   or \
#     matched(line,"^ *end program")  :
#     read=False

#def get_args(line):
#  global args,inarg
#  if inarg:
#    for arg in line.replace("subroutine","")          \
#                    .replace("function","")          \
#                    .replace(path.splitext(fich)[0],"")          \
#                    .replace( " ","") \
#                    .replace( "::","") \
#                    .replace( "&","") \
#                    .replace( "(","") \
#                    .replace( ")","")   \
#                    .replace("\n","").split(","):
#      if arg:
#        args.append(arg)

#def get_alloc(line):
#  global inalloc,allocated
#  if inalloc:
#    listline1=line.replace(" allocate(","") \
#                  .replace(" ALLOCATE(","")   \
#                  .replace("\n","")   \
#                  .replace(" ","").split("),")
#    for var in listline1:
#       allocated.append(var.split("(")[0])

#def get_var(line):
#  global indim,isdim,dim,listvar
#  listline1=line.replace("integer","")          \
#                .replace("double precision","") \
#                .replace("logical","")          \
#                .replace("allocatable","")          \
#                .replace(" dimension ","")      .replace( " ","") \
#                                                .replace( "::","") \
#                                                .replace( "&","") \
#                                                .replace("\n","").split(",")
#  indim=False
#  isdim=False
#  dim=[]
#  for var in listline1:
#    if var:
#      if "(" in var:
#        indim=True
#        var1=var.split("(")[0]
#        dim=[var.split("(")[1].replace(")","")]
#        if "dimension" in var1:
#          isdim=True
#      elif indim:
#        dim.append(var.replace(")",""))
#      if ")" in var:
#        if not "dimension" in var1:
#          if typepvar :
#            listvar[var1]=[typepvar,dim]
#          else:
#            print fich,var
#        indim=False
#        if not isdim:
#          dim=[]
#      elif not indim:
#        listvar[var]=[typepvar,dim]



##excluded_dim=["mtt","nmx","lgcmdx","ndimub","ndimctf","ndimnts","mdimub","mdimtbf","mdimtnf","mdimtcf","ip00","ip11","ip12","ip13","ip21","ip31","ip40","ip41","ip42","ip43","ip44"]
#excluded_dim=["mtb","mtt","nobj","nmx","lgcmdx","ndimub","ndimctf","ndimnts","mdimub","mdimtbf","mdimtnf","mdimtcf","ip00","ip11","ip12","ip13","ip21","ip31","ip40","ip41","ip42","ip43","ip44"]
### prettyfication
#listfiles = [ f for f in listdir(".") if isfile(join(".",f)) ]
##for fich in ["solve.f90"]:
#for fich in listfiles:
#  fileExtension = path.splitext(fich)[1]
#  if fileExtension == ".f90" :
##      print fich
#      number_of_units=count_units(fich)
#      i=0
#      read=False
#      i=0
#      while i < number_of_units:
#        i=i+1
#        listvar={}
#        k=0
#        read=False
#        typepvar=None
#        dimvar=False
#        inarg=False
#        inalloc=False
#        args=[]
#        allocated=[]
#        for line in open(fich):
#          if not matched(line,"^ *!") :
#            check_start(line,i)
#            get_args(line)
#            if read and not matched(line,"^ *!") :
#              if "allocate" in line or " ALLOCATE" in line :
#                inalloc=True
#              if "integer" in line :
#                typepvar="integer"
#              if "double precision" in line :
#                typepvar="real"
#              if "logical" in line :
#                typepvar="logical"
#              if matched(line,"[ ,]dimension[ (]"):
#                dimvar=True
#              get_alloc(line)
#              if typepvar or dimvar:
#                get_var(line)
#              if not "&" in line:
#                typepvar=None
#                dimvar=None
#                inalloc=False
#                inarg=False
#        k=0
#        read=False
#        decl=False
#        for line in fileinput.input(fich,inplace=1):
#          check_start(line,i)
#          if read and not matched(line,"^ *!") :
#            if "integer" in line or \
#                    "double precision" in line or \
#                    " dimension " in line or \
#                    "logical" in line or decl :
#              if "&" in line:
#                decl=True
#              else:
#                decl=False
#            else:
##              print "",
#              print line,

#            if "implicit none" in line:
#              integer = get_in_list(listvar,"integer",False,True,args)
#              if integer:
#                pprint("  integer          :: ",integer)
#              real = get_in_list(listvar,"real",False,True,args)
#              if real:
#                pprint("  double precision :: ",real)
#              logical = get_in_list(listvar,"logical",False,True,args)
#              if logical:
#                pprint("  logical          :: ",logical)

#              integer = get_in_list(listvar,"integer",True,True,args)
#              if integer:
#                pprint("  integer         ,allocatable :: ",integer)
#              real = get_in_list(listvar,"real",True,True,args)
#              if real:
#                pprint("  double precision,allocatable :: ",real)
#              logical = get_in_list(listvar,"logical",True,True,args)
#              if logical:
#                pprint("  logical         ,allocatable :: ",logical)

#              liste=get_in_list(listvar,"integer",True,True,args)
#              liste.extend(get_in_list(listvar,"real",True,True,args))
#              liste.extend(get_in_list(listvar,"logical",True,True,args))
#              for var in liste:
#                if not var.split("(")[0] in args and \
#                   not var.split("(")[0] in allocated:
#                  print "!WARNING : AUTOMATIC ALLOCATION\n allocate("+var+")"
#          else:
#            print line,


#stream=fileinput.input("compil.log")
#list_todo=[]
#for line in stream:
#  if "A statement function is an obsolescent feature in Fortran 2008" in line:
#    fich=line.split("(")[0]
#    line_numbr=line.split("(")[1].split(")")[0]
#    fx=stream.readline()
#    list_todo.append([fich,line_numbr,fx])
#for fich,line_numbr,fx in list_todo:
#  for newline in fileinput.input(fich,inplace=1):
#    if matched(newline,"^ *end subroutine") :
#      print "contains"
#      print "function"+fx.split("=")[0]
#      print "implicit none" 
#      print "integer ::"+fx.split("(")[0]+","+fx.split("(")[1].split(")")[0]
#      print fx.split("(")[0]+"="+fx.split("=")[1],
#      print "end function"
#    if line_numbr == str(fileinput.filelineno()):
#      print ""
#    elif "integer" in newline:
#      print newline.replace(" ","").replace("::"+fx.split("(")[0].replace(" ","")+",","::").replace(","+fx.split("(")[0].replace(" ","")+"\n","\n").replace(","+fx.split("(")[0].replace(" ","")+",",","),
#    else:
#      print newline,



#listfiles = [ f for f in listdir(".") if isfile(join(".",f)) ]
#for fich in listfiles:
#  fileExtension = path.splitext(fich)[1]
#  if fileExtension == ".f90" :
#        line1=""
#        for line in fileinput.input(fich,inplace=1):
#          if not (matched(line,"function") and matched(line1,"contains") and matched(line2,"end function")):
#            print line1,
#          line2=line1
#          line1=line
#        with open(fich,'a') as f: f.write(line1)

#remark #15300: LOOP WAS VECTORIZED
#remark #15301: BLOCK WAS VECTORIZED
#remark #15301: OUTER LOOP WAS VECTORIZED
#remark #15301: PARTIAL LOOP WAS VECTORIZED
#remark #15301: PERMUTED LOOP WAS VECTORIZED
#remark #15301: REVERSED LOOP WAS VECTORIZED
#remark #17100: DISTRIBUTED LOOP WAS AUTO-PARALLELIZED
#remark #17109: LOOP WAS AUTO-PARALLELIZED
#remark #17101: parallel loop shared={ .2.7_2upper_.11 } private={ } firstprivate={ mtcx mfc } lastprivate={ } firstlastprivate={ } reduction={ }

#listfiles = [ f for f in listdir("build.optrpt/CMakeFiles/solver_air.dir/") if isfile(join("build.optrpt/CMakeFiles/solver_air.dir/",f)) ]
#for fich in listfiles:
#  filename,fileExtension = path.splitext(fich)
#  if fileExtension == ".optrpt" :
#    line_loop=""
#    list_vect={}
#    for line in open("build.optrpt/CMakeFiles/solver_air.dir/"+fich):
#      if "LOOP BEGIN" in line:
#        line_loop=int(line.split("(")[1].split(",")[0])
#      if " LOOP WAS AUTO-PARALLELIZED" in line:
#        list_vect[line_loop]=([],[],[],[],[],[])
#      if "remark \#17101" in line:
#        list_vect[line_loop][0]=line.split("{")[1].split("}")[0].split(" ")
#        list_vect[line_loop][1]=line.split("{")[2].split("}")[0].split(" ")
#        list_vect[line_loop][2]=line.split("{")[3].split("}")[0].split(" ")
#        list_vect[line_loop][3]=line.split("{")[4].split("}")[0].split(" ")
#        list_vect[line_loop][4]=line.split("{")[5].split("}")[0].split(" ")
#        list_vect[line_loop][5]=line.split("{")[6].split("}")[0].split(" ")
#        list_vect[line_loop][6]=line.split("{")[7].split("}")[0].split(" ")
#    if list_vect:
#      print list_vect
#      for line in fileinput.input(filename+".f90",inplace=0):
#        if fileinput.filelineno() in enumerate(list_vect):
#          print list_vect
#        print line,


#listfiles = [ f for f in listdir(".") if isfile(join(".",f)) ]
#for fich in listfiles:
#  fileExtension = path.splitext(fich)[1]
#  if fileExtension == ".f90" :
#        test=False
#        test1=False
#        i=0
#        j=0
#        for line in fileinput.input(fich,inplace=0):
#          if test and "::" in line:
#            i=fileinput.filelineno()
#            test1=True
#          if test and "parameter" in line:
#            i=fileinput.filelineno()
#            test1=True
#          if test1 and "end subroutine" in line:
#            j=fileinput.filelineno()
#            break
#          if test1 and "end function" in line:
#            j=fileinput.filelineno()
#            break
#          if test1 and "return" in line:
#            j=fileinput.filelineno()
#            break
#          if test1 and "contains" in line:
#            j=fileinput.filelineno()
#            break
#          if matched(line,"^ *contains"):
#            test=True
#        fileinput.close()
#        for line in fileinput.input(fich,inplace=1):
#          if fileinput.filelineno() == j :
#            print "!$OMP END MASTER"
#          print line,
#          if fileinput.filelineno() == i :
#            print "!$OMP MASTER"



listfiles = [ f for f in listdir(".") if isfile(join(".",f)) ]
for fich in listfiles:
  fileExtension = path.splitext(fich)[1]
  if fileExtension == ".f90" :
        for line in fileinput.input(fich,inplace=1):
          if not "!$OMP" in line:
            print line,



