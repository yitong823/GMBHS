# GMBHS (A Pipeline for Gene Modeling based on HomologousHomologous Search)
##简介
    对于非模式生物的基因组草图序列来说，同源搜索是一项非常重要的注释信息来源，而且准确性很高，一般同源搜索获得的高相似位置都是真实的编码蛋白的某个基因所在，然而它有一个非常重要的缺陷，就是无法获得精确的基因结构，包括上下游和高相似内部区域的外显子内含子间接识别位点等。故而，我们在同源搜索的基础之上，基于python编写一款软件，利用同源搜索的结果，提取相应的基因组序列，进而使用基因建模软件对其进行结构建模，以获取精确的基因结构模型，最后将结果以gff3格式输出，可直接用于基因组可视化工具，如，IGV或JBrowse。
    
 
