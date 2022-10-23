#! /usr/local/bin/python3.7
# -*- coding: utf-8 -*-
"""
@author: 焱翊曈
"""
#模块导入检查准备工作
try:import os
except:print("模块导入失败：os,程序终止！");exit()
try:import argparse
except:print("模块导入失败：argparse,程序终止！");exit()

#程序参数及设定
parser = argparse.ArgumentParser(description='猴痘病毒基因组序列比对工具')
parser.add_argument('-p','--protein',metavar='',help='已知蛋白质序列文件(format:xxx.fasta)')
parser.add_argument('-g','--genome',metavar='',help='基因组序列文件 (format:xxx_genome.fasta)')
parser.add_argument('-o','--outdir',metavar='',help='结果输出目录')
parser.add_argument('--identity',type=float,metavar='',help='结果过滤中的一致性比例',default = 90)
parser.add_argument('--distance',type=int,metavar='',help='结果整合时的距离阈值',default = 100)

args = parser.parse_args()

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

#日志记录文件
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

#接收以及解析参数
def handle_parameter():
    try:import sys
    except:print("fail to import：sys");exit()
    try:import getopt
    except:print("fail to import：getopt");exit()

    try:
        opts, args = getopt.getopt(sys.argv[1:], "-h:-p:-g:-o::-f", ["help","protein=","genome=","out_dir=","blast_gff_file=",])
    except getopt.GetoptError:
        print("\n参数使用错误！\n")
        sys.exit(2)
    log(opts,logfile,1)

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
    

#创建结果输出目录并返回创建目录的路径
def output_file(file):
    if file ==None:
        os.system("mkdir output") #默认在当前工作目录下创建输出目录output
        output_path = cwork_dir+'/output'
        return output_path
    else:
        os.system("mkdir -p "+file)
        return file

#创建本地blast数据库
def make_blastdb(genome_file):
    
    os.system("mkdir "+ output+"/blastdb")
    title = genome_file.split("/")[-1].split("/")[0].split(".")[0]
    out = output+"/blastdb/"+title
    input_type = 'fasta'
    dbtype = 'nucl'
    command = "makeblastdb -in %s -input_type %s -title %s -dbtype %s -out %s "%(genome_file,input_type,title,dbtype,out)
    os.system(command)
    log("make the blastdb successfully!",logfile,3)
    return title

#调用tblastn进行比对,返回比对文件的绝对路径
def tblastn(protein_file,db_name):
    db = output+"/blastdb/"+db_name
    out = output +"/"+ db_name+"_tblastn.outfmt6"
    command = "tblastn -query %s -db %s -out %s -outfmt 6 -evalue 1e-5 "%(protein_file,db,out)
    os.system(command)
    log("blast alignment has done!",logfile,3)
       
    return out

#对比对结果进行过滤和区间整合
def HSP_handle(file,gap_len,identity):
    f = open(file,"r")
    infos = f.readlines()
    ff = open("./temp.txt",'w')
    for i in infos:
        if int(float(i.split("\t")[2]))<identity:
            infos.remove(i)
    for j in infos:
        idname= "_".join(j.split("\t")[:2])
        start = j.split("\t")[8]
        end = j.split("\t")[9]
        score = j.strip().split("\t")[11]
        list1 = [idname,start,end,score]
        realinfo = "\t".join(list1)+"\n"
        ff.write(realinfo)
    
    f.close()
    ff.close()
    
    #对处理成bed格式的文件使用bedtools进行合并
    #正链合并操作
    zheng_strand = "cat ./temp.txt | awk -v OFS='\t' '{if($2<$3)print $0}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge -d "+str(gap_len)+" -c 4 -i - -o sum | awk '{print $0 \"+\" }' > ./zheng_strand.txt"
    os.system(zheng_strand)
    #负链合并操作
    fu_strand = "cat ./temp.txt | awk -v OFS='\t' '{if($2>$3)print $1,$3,$2,$4}' | sort -k1,1 -k2,2n -k3,3n | bedtools merge -d "+str(gap_len)+" -c 4 -i - -o sum | awk '{print $0 \"-\" }' > ./fu_strand.txt"
    os.system(fu_strand)
    #正链和负链文件结果合并
    command = "cat ./zheng_strand.txt ./fu_strand.txt | sort -k2,2 -k3,3n > ./merge.txt"
    os.system(command)
    fff = open("./merge.txt",'r')
    HSPs = fff.readlines()
    os.remove("./temp.txt")
    os.remove("./zheng_strand.txt")
    os.remove("./fu_strand.txt")
    os.remove("./merge.txt")
    
    
    return HSPs   

    

#将整合结果转为gff3格式

def get_gff3(hsp_list,outfile):
    f = open(outfile,'w')
    hlist =hsp_list
    for i in hlist:
        gene = i.split("_")[0]
        chro = i.split("\t")[0].split("_")[1]
        start = i.split("\t")[1]
        end = i.split("\t")[2]
        lastcol = i.split("\t")[3].strip()
        if lastcol.endswith("+"):
            strand = "+"
            score = lastcol.strip("+")
        else:
            strand = "-"
            score = lastcol.strip("-")
        zhushi = "NAME="+gene+";"
        innfo = [chro,"TBLASTN","gene",start,end,score,strand,".",zhushi]
        info_line = "\t".join(innfo)+"\n"
        f.write(info_line)
    
    f.close()
    log("changed outfmt6 to gff3 format ,Done!",logfile,3)
    
    return None
    
            
            
        
        
        

#调用blast_gff.py文件将比对文件outfmt6格式转换为gff3格式
'''def blastgff(gfffile,outfmt_file,result_file):
    command = "python "+gfffile+" -b "+outfmt_file+" > "+ result_file
    os.system(command)
    log("changed outfmt6 to gff3 format ,Done!",logfile,3)

    return None'''
    
#生成基因列表文件
def get_gene_list(input_file,outfile):
    f = open(input_file,"r")
    genelist = []
    infos = f.readlines()
    for i in infos:
        if i.startswith(">"):
            gene = i.replace(">", "")
            if gene not in genelist:
                genelist.append(gene)
            else:
                continue
        else:
            continue
    f.close()
    ff = open(outfile,"w")
    for j in genelist:
        ff.write(j)
    ff.close()
    log("get the genelist file, Done!",logfile,3)
    
    return None


#运行主程序
def main(protein_file,db_name):
    #进行tblastn比对，并将结果输出到结果目录
    blast_out = tblastn(protein_file,db_name)
    HSPs_list = HSP_handle(blast_out,args.distance,args.identity)
    gfffile = output+"/tblastn_result.gff"
    get_gff3(HSPs_list, gfffile)
    genelist_file = output+"/genelist.txt"
    get_gene_list(args.protein, genelist_file)
    
    
    
    return None


if __name__ =='__main__':
    
    output = output_file(args.outdir) #创建结果输出目录
    logfile = output+"/processing.log"  #创建运行日志输出文件
    handle_parameter() #进行参数解析
    check_blast("blastp")
    check_blast("tblastn")
    check_blast("makeblastdb") #检查blast相关程序
    check_blast("bedtools")
    cwork_dir = os.getcwd()
    db_name = make_blastdb(args.genome)
    main(args.protein,db_name)



