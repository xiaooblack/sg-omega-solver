# sg-omega-solver V-0.1
Vertical velocity diagnosis of Omega equation under semi-geostrophic framework
## 背景
1. 基于准地转理论的Omega方程对垂向速度进行诊断已经较为成熟，但准地转框架下诊断出的垂向速度仅含有地转拉伸强迫，对于非地转成分显著的亚中尺度运动所伴随的次级环流难以准确刻画。 
2. Thmoas,2008中使用半地转框架对次级环流进行了解析，说明半地转框架具有描述亚中尺度次级环流的能力。  
3. 半地转框架由于缺乏非地转的局地时间变化，难以描述次级环流向下强化的现象。（缺少MC williams提出的advective feedback）  
## 目标
1. 在半地转框架中加入非地转的局地时间变化项，推导Hoskins形式的半地转Omega方程
2. 对方程进行离散，利用松弛迭代法进行对垂向速度进行诊断
## 作者的话
希望这份文档能对读者有所帮助，目前仅有matlab代码，python的开发可能要随缘。  
如果这份文档对你有帮助，请保佑作者原神抽卡不歪，万分感谢！！！！
