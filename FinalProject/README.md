12212706 刘一进

Final-Project 运行方法

所有代码相关文件均存放在 FinalProject 文件夹中
其中 main.m 为主程序
	 hole.msh 为带孔的网格
	 square.msh 为正方形不带孔的网格

默认上传到github的版本中，main.m 是对应于第三问manufactured problem 的代码
直接运行会输出误差值以及3张绘制的图像


设置 平面应力/应变 问题
	设置 problem_type 变量
		0为平面应力
		1为平面应变

设置source term
	analy_u 为u的解析解
	analy_du 为u的导数的解析解

设置应力边界条件(Neumann)
	sigma 为边界上应力对应函数
		sigma(1) 为xx应力
		sigma(2) 为yy应力
		sigma(3) 为xy应力

设置位移边界条件(Dirichlet)
	"g = zeros(n_np, 2);"
	通过给g赋值设定对应节点上的位移
		g(:, 1) 为节点x位移
		g(:, 2) 为节点y位移		

设置边界条件

	对不同的网格需要调整 getBoundary 函数
	该函数通过判断x,y坐标值，来划分边界
	对不同的网格需要调整其中的内容
	具体见代码注释

	对于同一类网格在76-99行处的switch中设置边界
		node_type = 0 内部节点
				  = 1 Dirichlet边界节点
				  = 2 Neumann边界节点
		ID(:, 1)控制x方向自由度
		ID(:, 2)控制y方向自由度
		例：
			x,y 都无约束
			"counter = counter + 2;
            ID(ii,1) = counter - 1;
            ID(ii,2) = counter;"

            x无约束，y有约束(位移值在 h 中设置)
			"counter = counter + 1;
            ID(ii,1) = counter;"

            x,y 都无约束
			不需要代码设置