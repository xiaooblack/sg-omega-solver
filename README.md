# sg-omega-solver
Vertical velocity diagnosis of Omega equation under semi-geostrophic framework
## 背景
基于准地转理论的Omega方程对垂向速度进行诊断已经较为成熟，但准地转框架下诊断出的垂向速度仅含有地转强迫，因此无法准确刻画亚中尺度次级环流。  
Thmoas,2008中使用半地转框架对次级环流进行了解析，说明半地转框架具有描述亚中尺度次级环流的能力。  
然而半地转框架由于缺乏非地转的局地时间变化，难以描述次级环流向下强化的现象。（缺少MC williams提出的advective feedback）  
## 目标
1. 在半地转框架中重新加入非地转的局地时间变化项，依照广义Omega方程重新进行Omega方程的推导
2. 对方程进行离散，利用松弛迭代法进行对垂向速度进行诊断
