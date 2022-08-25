function func = getFunc(C)
%
% Parametrized_benchmark_func evaluates (pop) Population on function (fun_n) using CEC2021 benchmark
% According to  the Parametrized Selection Vector (C), where:
% C1 is Bias indicator (0=0 ,1=F*)
% C2 Shift indicator (0=0, 1=oi)
% C3 Rotation indicator (0=I 1=M)



switch num2str(C)
    
    case '0  0  0'
        func=@cec21_basic_func;
    case '1  0  0'
        func=@cec21_bias_func;
    case '0  1  0'
        func=@cec21_shift_func;
    case '0  0  1'
        func=@cec21_rot_func;
    case '1  1  0'
        func=@cec21_bias_shift_func;
    case '1  0  1'
        func=@cec21_bias_rot_func;
    case '0  1  1'
        func=@cec21_shift_rot_func;
    case '1  1  1'
        func=@cec21_bias_shift_rot_func;
    otherwise
        disp('Undefined Selection Vector')
        
end
end