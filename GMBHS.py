#! /usr/local/bin/python3.7
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 20 19:57:46 2020

@author: 焱翊曈
"""

try:from Bio import SeqIO
except:print("模块导入失败：SeqIO(Bio),程序终止！");exit()
try:import os
except:print("模块导入失败：os,程序终止！");exit()
try:import argparse
except:print("模块导入失败：argparse,程序终止！");exit()
parser = argparse.ArgumentParser(description='同源搜索和基因结构建模批处理')
parser.add_argument('-p','--protein',metavar='',help='已知蛋白质序列文件(format:xxx.fasta)')
group = parser.add_mutually_exclusive_group()
group.add_argument('-g','--genome',metavar='',help='基因组序列文件 (format:xxx_genome.fasta)')
group.add_argument('-w','--genome_remote',metavar='',help=' 下载基因组序列的网址链接')
parser.add_argument('-o','--outdir',metavar='',help='结果输出目录')
parser.add_argument('-s','--step',type=int,metavar='',help='迭代提取基因组序列的步长',default= 50)
parser.add_argument('--score',type=int,metavar='',help='结果过滤中的相似性分值',default = 100)
parser.add_argument('--identity',type=float,metavar='',help='结果过滤中的一致性比例',default = 50)
parser.add_argument('--distance',type=int,metavar='',help='结果整合时的距离阈值',default = 1000)
parser.add_argument('--local_modeling_tool',metavar='',help='建模所使用的工具(默认为augustus)',default = 'augustus')
parser.add_argument('--max_iteration',type = int,metavar='',help='延伸序列最大迭代次数',default = 100)
parser.add_argument('--num_threads',type = int,metavar='',help='blast比对时设置的所用线程数',default = 1)
args = parser.parse_args()
    
def import_module(module_name, logfile = None):
    # 导入模块
    try:
        module = __import__(module_name) 
    except:
        log("模块导入失败:%s,程序终止！"%module_name,file = logfile)
        exit()
    else:
        return module



    # 执行系统命令，解析进程
def subprocess_popen(statement):
    try:import subprocess
    except:print("模块导入失败：subprocess,程序终止！");exit()
    p = subprocess.Popen(statement, shell=True, stdout=subprocess.PIPE) 
    while p.poll() is None:  
        if p.wait() != 0:  
            # print("命令执行失败:%s,程序终止！"%statement)
            return (False,p.returncode)
        else:
            re = p.stdout.readlines()  # 获取原始执行结果
            result = []
            for i in range(len(re)):  # 由于原始结果需要转换编码，所以循环转为utf8编码并且去除\n换行
                res = re[i].decode('utf-8').strip('\r\n')
                result.append(res)
            return (result,p.returncode)

#日志记录
def log(info,file, flag = 3):
    # 输出提示信息（info）
    # 1：输出到日志（file）
    # 2：输出到屏幕
    # 3：输出到日志和屏幕
    info_list = []
    if type(info) == str:
        info_list = [info]
    elif type(info) == list:
        for i in info:
            if type(i) == tuple:
                listinfo = list(i)
                prominfo = ":".join(listinfo)
                info_list.append(prominfo)
                info = ",".join(info_list)
            else:
                info_list.append(i.rstrip("\n"))
                info = "\n".join(info_list)

    
        
    if not file or flag == 2:
        print(info)
    elif flag == 1:
        with open(file,"a") as flog:
            flog.write(info + "\n")
    elif flag == 3:
        with open(file,"a") as flog:
            flog.write(info + "\n")
        print(info)

#获取当前工作路径
cwork_dir = os.getcwd()

#蛋白质fasta格式数据校验和读取
def check_get_protein(protein_file):
    if not protein_file.endswith(".fasta" or ".fa"):
        log("format error: the file must be the format of '.fasta' or '.fa'! \n Processing termination!",logfile,1)
        exit()
    elif protein_file.endswith(".fasta" or ".fa" ):
        file = cwork_dir+"/"+protein_file
        if os.path.exists(file) ==True:
            print("the file is exist!")
            pro_list = []
            for info in SeqIO.parse(file,"fasta"):
                protein_id = str(info.id)
                pid = protein_id.split("|")[1]
                protein_seq = str(info.seq)
                information = ">"+pid+"\n"+protein_seq
                pro_list.append(information)
                log("check the "+ pid+" sequence successfully!",logfile,2) 
            return pro_list
        else:
            print("not exist the file!")
    else:
        log("No found protein file \n Processing termination!",logfile)
        exit()


#基因组序列文件格式校验
def check_genomefile(genome_file):
    # 返回0.则检验有误
    #返回1，则检验成功
    if not genome_file.endswith(".fasta" or ".fa"or ".fna"):
        log("the format of genome file is error!\n please reupload the file or use '-g' to get file online ",logfile)
        return 0
    else:
        log("check the genome file successfully!",logfile)
        return 1
#文件夹校验
def dir_check(dir):
        try:from pathlib import Path
        except:print("模块导入失败：pathlib.Path");exit()
        return Path(dir).is_dir()

def file_type(file):
    #返回文件类型
    return file.rsplit(".",1)[-1]

def file_check(file, logfile = None, ftype = None):
    #检查文件是否合法
    # -1:文件不存在；0:文件不合法；1:文件存在、合法
    os = import_module("os", logfile)
    if ftype:
        if os.path.isfile(file):
            if type(ftype) == str:
                if file_type(file) == ftype:
                    return 1
                else:
                    return 0
            elif type(ftype) == tuple:
                if file_type(file) in ftype:
                    return 1
                else:
                    return 0 
        else:
            return -1
    else:
        if os.path.isfile(file):
            return 1
        else:
            return -1

def com_file_check(file, dc_type, logfile = None):
    #判断是否为某个类型的压缩文件
    re = import_module("re", logfile)
    # file = file.rsplit("/",1)[-1]
    obj = re.match("(.*)\\" + dc_type,file)
    if not obj:
        return False
    else:
        return obj.group(1)

def py_dc(file,dc_type, outdir, logfile = None):
    # 使用python的gzip、tarfile、zipfile库解压
    # 1:解压成功；0:不支持dc_type压缩格式；e:解压错误
    os = import_module("os", logfile)
    try:
        if dc_type == ".tar.gz":
            gzip = import_module("gzip", logfile)
            tarfile = import_module("tarfile", logfile)
            
            tar_file = file.rstrip(".gz")
            gz_file = gzip.GzipFile(file)
            with open(tar_file, "wb+") as ft:
                ft.write(gz_file.read())
            gz_file.close()
            
            tar = tarfile.open(tar_file)
            tar.extractall(outdir)
            tar.close()
            os.remove(file)
            os.remove(tar_file)
            return 1
        elif dc_type == ".gz":
            gzip = import_module("gzip", logfile)
            
            file_name = file.rstrip(".gz")
            gz_file = gzip.GzipFile(file)
            with open(file_name, "wb+") as ft:
                ft.write(gz_file.read())
            gz_file.close()
            os.remove(file)
            return 1
        elif dc_type == ".tar":
            tarfile = import_module("tarfile", logfile)

            tar = tarfile.open(file)
            tar.extractall(outdir)
            tar.close()
            os.remove(file)
            return 1
        elif dc_type == ".zip":
            zipfile = import_module("zipfile", logfile)

            zip_file = zipfile.ZipFile(file)
            zip_file.extractall(outdir)
            zip_file.close()
            os.remove(file)
            return 1
        elif dc_type == ".rar":
            try:
                from unrar import rarfile 
            except:
                log("模块导入失败:from unrar import rarfile,程序终止！", logfile)
                raise Exception("")
          
            rar_file = rarfile.RarFile(file)
            rar_file.extractall(outdir)
            os.remove(file)
            return 1
        else:
            return 0
    except Exception as e:
        return e


def decompression(file, os, logfile = None):
    #解压缩文件
    types = [".tar.gz",".tar.xz",".tar.bz2",".zip",".rar",".7z",".tar",".gz",".bz2",".Z"]    
    commands = ["tar zxvf","tar -xvJf","tar jxvf","unzip","rar x","7z x","tar xvf","gzip -d","bzip2 -d","uncompress"]
    if file_check(file) == 1:      # file存在时
        for i in range(len(types)):
            dc_file = com_file_check(file,types[i])
            if dc_file and file_check(dc_file) == 1:  # 不能存在已经解压好的文件,否则解压会失败
                log("已经存在解压好的文件:%s"%dc_file, logfile)
                return dc_file
            if dc_file:  
                if os == "poxsix":
                    command = commands[i] + " " + file
                    result,returncode = subprocess_popen(command)
                    if not returncode:
                        log("基因组文件解压完成:%s\n"%file, logfile)
                        return dc_file
                    # elif returncode == 127:
                    else:
                        log("使用系统命令%s解压失败:"%command, logfile)
                        log(result, logfile)
                        log("使用python库函数解压.", logfile)
                        dc_check = py_dc(file,types[i],output)
                        if dc_check == 1:
                            log("基因组文件解压完成:%s\n"%dc_file, logfile)
                            return dc_file
                        elif dc_check == 0:
                            log("不支持该压缩格式:%s,程序终止！"%file, logfile)
                            raise Exception("")     
                        else:
                            log("解压失败:%s,程序终止！"%dc_check, logfile)
                            raise Exception("")
                else:
                    dc_check = py_dc(file,types[i],output)
                    if dc_check == 1:
                        log("基因组文件解压完成:%s\n"%dc_file, logfile)
                        return dc_file
                    elif dc_check == 0:
                        log("不支持该压缩格式:%s,程序终止！"%file, logfile)
                        raise Exception("")
                    else:
                        log("解压失败:%s,程序终止！"%dc_check, logfile)
                        raise Exception("")                      
    else:
        log("压缩文件不存在:%s,程序终止！"%file, logfile)
        raise Exception("")

#参考基因组数据库校验
def check_blastdb(db):
    if db[-1] == "/": #不能是一个目录
        return 0
    else:
        db_dir = db.rsplit("/",1)[0]
        db_title = db.rsplit("/",1)[1]
        if dir_check(db_dir):
            files,returncode = subprocess_popen("ls " + db_dir)
            flag = 0
            for file in files:
                if file.split(".")[0] == db_title and file.split(".")[1] in("nin","nhr","nsq"):
                    flag+=1
            if flag >= 3:
                log("db check successfully!",logfile,1)
                return db_title
            else:
                log("not have all of 'nin','nhr','nsq'",logfile,1)
                return db_title
        else:
            log("db check error!the db_format has problems",1)
            return db_title


 #blast相关程序校验

def check_blast(program):
    result,returncode = subprocess_popen(program + " -version")
    if returncode == 127:
        log("未安装%s或环境变量未配置。如果是Ubuntu系统，输入‘sudo apt install ncbi-blast+’安装。"%program,logfile)
        exit()
    elif returncode:
        log("%s无法使用,请检查权限。若没有访问权限，请联系管理员"%program,logfile)
        exit()
    else:
        log(result[0],logfile)
    
#samtools程序校验
def check_samtools(program):
    result,returncode = subprocess_popen(program + " --version")
    if returncode == 127:
        log("未安装%s或环境变量未配置。如果是Ubuntu系统，输入‘sudo apt install samtools’安装。"%program,logfile)
        exit()
    elif returncode:
        log("%s无法使用,请检查权限。若没有访问权限，请联系管理员"%program,logfile)
        exit()
    else:
        log(result[0],logfile)
       

#创建本地blast数据库
def make_blastdb(genome_file):
    
    os.system("mkdir "+ output+"/blastdb")
    title = genome_file.split("/")[-1].split("/")[0].split(".")[0]
    out = output+"/blastdb/"+title
    input_type = 'fasta'
    dbtype = 'nucl'
    command = "makeblastdb -in %s -input_type %s -title %s -dbtype %s -out %s "%(genome_file,input_type,title,dbtype,out)
    os.system(command)
    log("make the blastdb successfully!",logfile,2)
    return title


#创建结果输出目录并返回创建目录的路径
def output_file(file):
    if file ==None:
        os.system("mkdir output") #默认在当前工作目录下创建输出目录output
        output_path = cwork_dir+'/output'
        return output_path
    else:
        os.system("mkdir "+file)
        return file


#调用tblastn进行比对
def tblastn(protein_file,db_name,num_threads,genome_db):
    info = SeqIO.read(protein_file, "fasta")
    protein_name = str(info.id)
    outfile = os.listdir(output)
    if "blastdb" in outfile:
        db = output+"/blastdb/"+db_name
        out = output +"/"+ db_name +"-"+protein_name+ "_tblastn.outfmt6"
        command = "tblastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 -num_threads %s"%(protein_file,db,out,num_threads)
        os.system(command)
    else:
        db = genome_db
        out = output +"/"+ db_name +"-"+protein_name+ "_tblastn.outfmt6"
        command = "tblastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 -num_threads %s"%(protein_file,db,out,num_threads)
        os.system(command)
    return out

#对tblastn的结果给出建议
def blast_suggest(outfmt6_file):
    f = open(outfmt6_file,"r")
    HSPs = f.readlines()
    for i in HSPs:
        info = i.split("\t")
        if float(info[2])<20:
            log(info[0]+"_"+info[1]+"_"+info[6]+"-"+info[7]+":low indentity warning!",logfile)
        if float(info[11])<100:
            log(info[0]+"_"+info[1]+"_"+info[6]+"-"+info[7]+":low score warning!",logfile)


#多个 HSPs 的情况《=基因家族 情况对format6格式解析
def HSPs_handle(file,gap_len = 1000,identify = 50,score = 100):
    chr = {}
    HSP = []
    with open(file) as f:
        for line in f:
            line = line.split("\t")
            if float(line[2]) > identify or float(line[11]) > score:
                if line[1] not in chr:
                    chr[line[1]] = [[],[]]
                if line[8] < line[9]:
                    chr[line[1]][0].append([line[1],int(line[8]),int(line[9]),float(line[11])])
                else:
                    chr[line[1]][1].append([line[1],int(line[8]),int(line[9]),float(line[11])])
    for c in chr:
        for i in (0,1):
            if chr[c][i]:
                # print(chr[c][0])
                chr[c][i] = sorted(chr[c][i], key = lambda x : x[i+1])   # i=0,即正链的话，起始位置在第二为即i+1
                # print(chr[c][0])
                start = chr[c][i][0][i+1]     #0-1  1-2
                end = chr[c][i][0][2-i]         #0-2   1-1
                allscore = chr[c][i][0][3]
                if not chr[c][i][1:]:
                    # HSP.append([c, start, end, allscore])       # 正链是start-end，负链为end-start
                    HSP.append([c] + [start, end][::1-2*i] + [allscore])  # y=1-2x  使得i为0时是正序切片，i为1时为倒序切片
                else:
                    for hsp in chr[c][i][1:]:
                        if hsp[1] - end > gap_len:
                            HSP.append([c] + [start, end][::1-2*i] + [allscore])
                            start = hsp[i+1]
                            allscore = hsp[3]
                        else:
                            allscore += hsp[3]
                        end = hsp[2-i]
                    HSP.append([c] + [start, end][::1-2*i] + [allscore])  # 考虑最后的一个HSP

    return HSP


#从基因组中提取序列【设置迭代步长参数 step】并返回得到的fasta文件路径
def extract(proteinname,genome_file,chro,start,end,step):
    fasta = genome_file+" "
    #正链情况
    if start<end:
        final_start = start - step
        final_end = end + step        
        command = "samtools faidx "+ fasta + chro + ":" + str(final_start) + "-" + str(final_end)
        subseq = os.popen(command)
        info = subseq.readlines()
        name = chro+"_"+str(final_start) + "-" + str(final_end) 
    #负链情况    
    if start > end:
        final_start = end -step
        final_end = start +step
        command = "samtools faidx "+ fasta + chro + ":" + str(final_start) + "-" + str(final_end)          
        subseq = os.popen(command)
        info = subseq.readlines()
        name = chro+"_"+str(final_end) + "-" + str(final_start)
    subseq_path = output+'/'+proteinname+"_HSP/"+name+".fasta"
    f = open(subseq_path,"w")
    f.writelines(info)
    f.close()
    return subseq_path
    
#对建模工具进行检查
def check_modeling_tool(program):
            result,returncode = subprocess_popen(program + " --version")
            if returncode == 0:
                log("augustus checked successfully!",logfile)
            else:
                log("未安装%s或环境变量未配置。请从官网下载后使用(URL=http://augustus.gobics.de/)"%program,logfile)
                exit()
            
                
#使用Augustus进行基因结构建模 返回gff3文件路径
def augustus(pname,fasta,output):
    name = fasta.split('/')[-1].split('_')[0]
    outfile = output+"/"+pname+"_"+name +"_augustus_out.gff3"
    command = "augustus --gff3=on --outfile="+outfile+ " --species=saccharomyces_cerevisiae_S288C "+fasta
    os.system(command) 
    return outfile

#在Augustus建模结果的gff3文件中记录相关信息（预测基因命名，score值记录，记录蛋白质id，迭代次数记录）以及正负链去冗余
def gff_plu_score(gfffile,fasta,tblast_score,proteinid,times):
    info = fasta.split("/")[-1].split("_")[1].split(".")[0].split("-")
    chrome = fasta.split("/")[-1].split("_")[0]
    head = int(info[0])
    tail = int(info[1])
    f = open(gfffile,"r")
    all_info = f.read()
    infolist = all_info.split("###\n")   #以###分割建模的基因信息
    ff = open("temp1.gff3","w+")
    if head - tail < 0: #正链信息
        for i in infolist:
            if "\t+\t" in i:
                ff.write(i)
                ff.write("###\n")
            elif "\t-\t" in i:
                continue
            else:
                ff.write(i)
    else:  #负链信息
        for i in infolist:
            if "\t-\t" in i:
                ff.write(i)
                ff.write("###\n")
            elif "\t+\t" in i:
                continue
            else:
                ff.write(i)
        
    f.close()
    ff.close()
    #注释行添加比对分值,protein_id,迭代次数等信息
    ff = open("temp1.gff3","r")
    fff = open("temp2.gff3","w")
    data = ff.readlines()
    for i in data:
        if not i.startswith("#"):
            if "gene" in i:
                geneinfo = i.strip().split("\t")
                augustu_head = int(geneinfo[3])
                if head<tail: #正链结果
                    fact_start = head+augustu_head-1  #建模结果实际的起始位置：原始区间+augustus建模区间-1          
                    gene_id = "XLOC_"+chrome+"_"+str(fact_start)
                    a = i.strip()[:-2]
                    final_info = a+gene_id+";score="+str(tblast_score)+";pre_protein="+proteinid+";#times="+str(times)+"\n"
                    fff.write(final_info)
                else:
                    fact_start = tail+augustu_head-1  #建模结果实际的起始位置：原始区间+augustus建模区间-1          
                    gene_id = "XLOC_"+chrome+"_"+str(fact_start)
                    a = i.strip()[:-2]
                    final_info = a+gene_id+";score="+str(tblast_score)+";pre_protein="+proteinid+";#times="+str(times)+"\n"
                    fff.write(final_info)
            else:
                fff.write(i)
        else:
            fff.write(i)
    ff.close()
    fff.close()
    
    os.remove("temp1.gff3")
    os.remove(gfffile)
    os.rename("temp2.gff3",gfffile)
  
    return gfffile


#解析建模结果：
def augu_gff(file):
    f = open(file,"r")
    lines = f.readlines()
    info = []
    for i in lines:
        if not i.startswith("#"):
            line = i.split("\t")
            info.append(line[2])
    if "start_codon" in info  and "stop_codon" in info:
        return True
    else:
        return False
#从augustus预测结果中提取蛋白质序列
def pred_protein(gff_file,output):
    f = open(gff_file)
    info = f.readlines()
    name = gff_file.split("/")[-1].split('augustus')[0]
    protein_list = []
    for i in info:
        i.strip(" ")
        if i.startswith("# protein sequence"):
            diyihang = i.split('[')[1].strip()
            if ']' in diyihang:
                protein = diyihang.split(']')[0]
                protein_list.append(protein)
                break
            else:     
                protein_list.append(diyihang)
                hang = info.index(i)
                while True:
                    line = info[hang+1].strip()
                    if  not line.endswith("]"):
                        hang_seq = line.split(" ")[1]
                        protein_list.append(hang_seq)
                        hang+=1
                        continue
                    else:
                        hang_seq = info[hang+1].split(" ")[1]
                        hang_seq = hang_seq.strip()
                        hang_seq = hang_seq.strip(']')
                        protein_list.append(hang_seq)
                        break
    preid = ">"+name+"\n"                  
    protein_list.insert(0, preid)
    protein_seq = "".join(protein_list)+"\n"
    f.close()
    protein_file = open(output+'/'+name+"_pred.fasta",'w')
    protein_file.write(protein_seq)
    protein_file.close()
    out = output+"/"+name+"_pred.fasta"
    return out

#使用blastp进行protein双序列对比：
def blastp(subject_file,protein_file):
    query = protein_file
    subject = subject_file
    name = subject.split("/")[-1].split('_pred')[0]
    out = output+"/"+ name + "_blastp.outfmt6"
    command = "blastp -query %s -subject %s -out %s -outfmt 6"%(query,subject,out)
    os.system(command)
    return out
#提取blastp的比对结果分值
def get_score(file):
    f = open(file,"r")
    info = f.readlines()
    allscore = 0
    for i in info:
        score = float(i.split("\t")[-1])
        allscore+=score
    return allscore
    
def remote_db(dir):
    dirs=os.listdir(dir)
    for file in dirs:
        if file.endswith(".fasta" or ".fa"):
            return file

#接收以及解析参数
def handle_parameter():
    try:import sys
    except:print("fail to import：sys");exit()
    try:import getopt
    except:print("fail to import：getopt");exit()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "-h-p:-db:-g:-w:-o:-s:", ["help","protein=","genome_db=","genome=","genome_remote=","out_dir=", \
                                                       "step=","score=","identity=","distance=","local_modeling_tool=","max_iteration=",])
    except getopt.GetoptError:
        print("\n参数使用错误！\n")
        sys.exit(2)
    log(opts,logfile,1)

def genomefile(cwork_dir):
    #两参数互斥，若没有创建好本地blast库，便使用基因组数据进行创建
    os.system(" cd "+cwork_dir)

    if args.genome:
        check_result=check_genomefile(args.genome)
        if check_result == 1:
            db_name = make_blastdb(args.genome)
            return db_name,args.genome
        else:
            exit()
            log("检查genome文件失败有误，退出程序！",logfile)

    elif args.genome_remote:
        os.system("wget -P "+output+" "+args.genome_remote)
        genome_file = os.path.join(output,os.path.basename(args.genome_remote))
        log("下载基因组序列文件中...",logfile)
        log("开始下载基因组序列文件。",logfile,1)
        genome_file_check = file_check(genome_file, logfile,("fasta", "fna", "fa"))
        if genome_file_check == -1:
            log("下载的基因组序列文件不存在:%s,程序终止！"%genome_file,logfile)
            exit()
        elif genome_file_check == 0:
            #文件解压
            log("基因组序列文件已下载:%s"%genome_file,logfile)
            download_genome = decompression(genome_file, os.name)                  
            if download_genome:
                if file_check(download_genome, logfile,("fasta", "fna", "fa")) == 0:
                    log("基因组序列文件不是fasta格式:%s,程序终止！"%download_genome,logfile)
                else:
                    log("基因组序列解压和检查成功",logfile)
                    db_name = make_blastdb(download_genome)
                    return db_name,download_genome
            else:
                log("不支持该压缩格式:%s,程序终止！"%genome_file,logfile)
                exit()
        elif genome_file_check == 1:
            db_name = make_blastdb(genome_file)
            log("基因组序列检查成功！",logfile)
            return db_name, download_genome


def main(outpath,protein,db_name,genome_file):

    result = check_get_protein(protein)
    gff3 = open(outpath+"/"+db_name+"_tblastn_augustus.gff3","a+") 
    pred_proteins_file = open(outpath+"/"+db_name+"_tblastn_augustus_pred.fasta","a+")
    all_blastp = open(outpath+"/"+db_name+"_tblastn_augustus_blastp.outfmt6","a+")

    for info in result:
        f = open("tem_protein.fasta","w")
        f.write(info)
        f.close()


        protein_file = os.path.realpath("tem_protein.fasta")

        db = db_name
        num_threads = args.num_threads #获取blast比对的线程数

        '''file_list = os.listdir()
        for i in file_list:
            if i.endswith("genome.fasta"):
                break
        genome_file = cwork_dir+'/'+i #获取基因组序列文件'''

        #进行tblastn比对，并将结果输出到结果目录 
        blast_out = tblastn(protein_file,db,num_threads,db_name) 
         
        blast_suggest(blast_out) #对blast结果给出建议
        HSP_list = HSPs_handle(blast_out,args.distance,args.identity,args.score) #解析比对结果，并提取独立的HSPs
        new_HSP_list = []
        for yyt in HSP_list:
            new_info = [yyt[0].replace("_","."),yyt[1],yyt[2],yyt[3]]
            new_HSP_list.append(new_info)

        proname = info.split('\n')[0].split(">")[-1]
        HSP_file = output+"/"+proname+'_'+"HSP"
        os.system('mkdir '+HSP_file) #创建output下的保存提取序列的子文件夹
        for HSP in new_HSP_list:
            subseq = extract(proname,genome_file,HSP[0],HSP[1],HSP[2],0) #提取序列并得到fasta路径
            gff = augustus(proname,subseq, output) #使用Augustus进行基因结构建模 *out.gff3
            if augu_gff(gff) ==True:  #解析建模结果是否完整
                continue
            else:
                step = 0
                flag = 1
                while augu_gff(gff) == False: #当建模结果不完整
                    if flag < args.max_iteration: #设置最大迭代次数

                        if min(int(HSP[1]),int(HSP[2]))-step-args.step<0:
                            break
                        else:
                            os.system("rm "+gff)
                            subseq = extract(proname,genome_file,HSP[0],HSP[1],HSP[2],step+args.step)
                            gff = augustus(proname,subseq, output) 
                            step+=args.step  #迭代提取序列
                            flag+=1
                    else:
                        log("超出最大迭代次数，系统程序终止或进入下一个蛋白建模！",logfile)
                        break

            #把每一个HSP结果追加到总的文档中  

            gff = gff_plu_score(gff,subseq,HSP[3],proname,flag) #在gff文档中记录blast比对相似性分值

            pred_protein_seq = pred_protein(gff, output) #提取Augustus预测到蛋白质序列,*pro.fasta

            


            blastpout = blastp(pred_protein_seq,protein_file) #进行blastp比对 *blastp.outfmt6

            fblastp = open(blastpout)
            blastp_info = fblastp.read()
            all_blastp.write(blastp_info)
            fblastp.close()

            final_name = blastpout.split("/")[-1].split('_blastp')[0]
            blastp_score = get_score(blastpout) #提取blast比对后的打分结果，并与同源搜索记录的分值进行对比
            if blastp_score >= float(HSP[3])*1: #阈值系数（暂定为1）
                log(final_name+":达到阈值要求，建模成功",logfile)
                HSP_gff = open(gff)
                HSP_gff_info = HSP_gff.read()
                gff3.write(HSP_gff_info)
                HSP_gff.close()

                pre_protein = open(pred_protein_seq)
                pre_protein_info = pre_protein.read()
                pred_proteins_file.write(pre_protein_info)
                pre_protein.close()
            else:
                log(final_name+":未达到阈值要求，建模结果warning！",logfile)
            os.system("rm "+ gff)
            os.system("rm "+ pred_protein_seq)
            os.system("rm " + blastpout)

            continue

        os.system('rm ' + protein_file)

    gff3.close()
    pred_proteins_file.close()
    all_blastp.close()

    


if __name__ =='__main__':
    
    output = output_file(args.outdir) #创建结果输出目录
    logfile = output+"/GMBHS.log"  #创建运行日志输出文件
    handle_parameter() #进行参数解析
    check_blast("blastp")
    check_blast("tblastn")
    check_blast("makeblastdb") #检查blast相关程序
    check_samtools("samtools") # 检查samtools程序
    check_modeling_tool(args.local_modeling_tool)#检查从头建模工具
    cwork_dir = os.getcwd()
    #获取当前目录下的文件名
    genomeinfo = genomefile(cwork_dir) #进行基因组检查，建库等
    db_name = genomeinfo[0]
    genome_file = genomeinfo[1]
    main(output,args.protein,db_name,genome_file)


