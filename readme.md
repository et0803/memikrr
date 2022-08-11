# 1. Aim
To replicate natural human arm posture by performing muscle effort minimization when given a desired hand(wrist) pose. The manuscript entitled ***Muscle-effort-minimization-inspired Kinematic Redundancy Resolution for Replicating Natural Posture of Human Arm*** has been accepted for publication in IEEE Transactions on Neural Systems and Rehabilitation Engineering.

<img src="https://github.com/et0803/memikrr/raw/main/doc/armSwivelAngleRange.gif" alt="armSwivelAngleRange" width="400">    <img src="https://github.com/et0803/memikrr/raw/main/doc/optimizedArmSwivelAngle.png" alt="optimizedArmSwivelAngle" width="400">

Figure note: a warmer color (red) in the muscles indicates a higher muscle effort (activation).

# 2. Dependencies
- Eigen 3.3
- Opensim 4.3
- orocos_kdl 1.5
- kdl_parser 1.13
- C++ 11

# 3. Test environment
- Ubuntu 18.04 (Linux 5.4.0)
- gcc version 7.5.0 (x86_64)
- cmake 3.18.4

# 4. Test commands
1. cd path_to_CMakeLists.txt
2. cmake . -DCMAKE_BUILD_TYPE=Release
3. make
4. ./memikrr ./data/Setup_Analyze.xml ./data/x_arm_musculoskeletal_lql.urdf ./data/x_arm_musculoskeletal_simplified_joints_aligned_lql_noJointOffset.urdf 2.0 0.01 0.1 0.1 2 ./data/wristPoseInput.txt ./result.txt
## 4.1 Argument list
1. relative file path of static optimization setup for arm osim model
2. relative file path of musculoskeletal arm urdf model
3. relative file path of simplified SRS arm urdf model 
4. optimizaiton initial guess interval of arm swivel angle in degree
5. local optimization stop criterion kappa
6. global optimization stop criterion1 epsilon
7. global optimization stop criterion2 delta
8. the index of the selected global stop criterion
9. relative file path of wrist poses
10. relative file path of optimization result output

# 5. Bug report
1. Each try of muscle effort estimation will cause a memory leakage of about 1M bytes  by calling API to Opensim AnalyzeTool StaticOptimization module.

# 6. Citation
Wait for production processing...