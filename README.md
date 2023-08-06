1 对kegg中所有与l1000中相关的drug的pathway进行整理
  先运行l1000_db脚本，再运行 cmap_db脚本，最后运行 kegg 脚本（pathway ground truth）

![gt](https://github.com/zxws1/insilicoMedicine_code/assets/12165605/75cd79d2-ad9f-40c8-99c9-1c3edd3b367e)


2 运行一个例子 Selumetinib

模型预测一个结果

![b](https://github.com/zxws1/insilicoMedicine_code/assets/12165605/401395b5-4041-4c8d-80b9-f52098c4f4e3)

对FDR的P-value进行升序排列，并与pathway ground truth中对应的cpd的pathway进行对比

![a](https://github.com/zxws1/insilicoMedicine_code/assets/12165605/733f0d6b-d2b5-46b7-95ae-bd61fa3cac22)
