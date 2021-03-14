# Genome-Annotation
同源基因指导的全基因组基因预测及结构建模

编写python脚本FilterGFF3.py，筛选augustus预测结果中不完整的基因。读入gff3结果文档，已gene为单位进行甄别操作：
	将预测正负链错误的gene删去；
	对未被预测出的序列删去，并在predicate_error.txt中进行标记；
	对被预测出多个匹配基因的序列，选择其与原序列overlap最大的一个保留，且若最大overlap < 0.6，同样删去，并标记；
	记录每一次迭代过程中未被完整预测的基因数目。

每次使用augustus软件进行基因预测的结果命名为augustus_out_*-0.gff3，将所有未被完整预测的基因信息另存为augustus_out_*-2.gff3文件，预测完整的数据另存为augustus_out_*-1.gff3文件（*表示迭代的次数，如此命名，方便进行批量迭代操作），并将每次预测成功的结果追加写入augustus_out_sum.gff3。

代码见附录1、FilterGFF3

对上述进⾏封装、迭代和预测
1、封装模块
1.1、接入“序列提取”模块
感谢『罗晓琦』同学提供的“提取序列”R语言脚本，现使用python调用

1.2、接入“判断预测结果完整性”模块
由于FilterGFF3.py已被封装完毕，只需导入模块并实例化对象，即可使用

