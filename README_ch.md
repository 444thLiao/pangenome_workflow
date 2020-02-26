# 泛基因组分析流程中文文档

本仓库是一个集成了多个软件，并通过luigi进行流程管理的综合性的流程。
概括性的讲，该流程接受的输入是pair end的二代测序单菌株全基因组数据，但也可以接受single end或者从网上下载基因组数据（但需要在data_input.tab中进行区分以上两种/多种数据类型)，通过完整的流程后输出多种结果。

其中包括
> 拼接后的基因组序列
> 注释后的蛋白质序列
> 结合拼接后的数据和原始数据推断的质粒的基因序列
> 拼接前后的质量、基因组基本信息的评估
> 传统的pubmed的mlst分型结果
> 基于16s/kraken2的物种注释
> 基于abricate的毒力因子/耐药因子的注释
> (Insertion sequence)IS区间的注释
> 基于roary的泛基因组分析结果

并基于以上的结果，导入geneious、plotly等可视化软件进行浏览以上的结果

总的来说，该流程执行的目的和使用的软件如下：

1. 质量评估(quality accessment) ([fastqc](https://github.com/s-andrews/FastQC))
2. 质量控制(quality control) ([trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic))
3. 综合性的质量评估([multiqc](https://github.com/ewels/MultiQC))
4. 拼接 ([shovill](https://github.com/tseemann/shovill)+spades)
5. 蛋白注释([prokka](https://github.com/tseemann/prokka#installation))
6. 物种注释([kraken2](https://github.com/DerrickWood/kraken2))
7. 基因组基本信息的评估([seqtk](https://github.com/lh3/seqtk))
8. 毒力因子/耐药因子的注释([abricate](https://github.com/tseemann/abricate))
9. (Insertion sequence)IS区间的注释([ISEscan](https://github.com/444thLiao/ISEScan/tree/test))
10. 基因组间物种距离的估计，以进行泛基因组分析前的预聚类(mash)
11. 泛基因组分析(roary)



#TODO list
1. 补充基于EZbio的16s物种注释，并入luigi流程中
2. 补充更详细的文档
3. 完善api部分，以实现分块完成
4. 加入更为丰富的流程定制化策略



