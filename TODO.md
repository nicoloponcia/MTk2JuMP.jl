
- say that custom lb,ub should be NaN in ScalesBounds.jl, then overwritten later. E.g. x,y for vehicle optimization
- add settings for control_der and control_acc specifying the type of finite difference
- better solver info handling. Rethink if possible to capture iptopt outputs without .txt
- recheck where use NonlinearExpr... and subsitute with general AbstractJuMPScalar