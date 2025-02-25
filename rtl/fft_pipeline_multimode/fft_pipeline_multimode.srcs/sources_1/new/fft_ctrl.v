`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: UCAS
// Engineer: Astron-fjh
// 
// Create Date: 2025/02/15 01:21:36
// Design Name: FFT Control Module
// Module Name: fft_ctrl
// Project Name: fft_pipeline_multimode
// Target Devices: 
// Tool Versions: 
// Description: 
//   This module controls the operation of the FFT (Fast Fourier Transform) process.
//   It manages the flow of data through the stages of the FFT computation by generating 
//   appropriate read and write addresses, enabling data read and write operations, and 
//   managing the state of the FFT process. It also controls the start and completion signals 
//   for the FFT computation. The module includes control for multiple stages of the FFT, 
//   including updating the stage number, group number, and butterfly operation parameters.
//   It interfaces with the SRAM and butterfly modules to complete the FFT calculations.
//
//   The module also generates control signals for the butterfly operation, including 
//   inputs for the real and imaginary parts of the input data, as well as the rotation 
//   factors for the FFT calculation. The FFT process is split into stages and operates 
//   on a set of data groups, with the control logic updating the stage and group counts 
//   accordingly. 
// 
// Dependencies:
//   - butterfly module for performing the butterfly operations in the FFT.
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module fft_ctrl(
// ============== 时钟信号 ===============
    input  wire               clk,
    input  wire               rst_n,

// ============== 输入信号 ===============
    input  wire               resort_complete,  // 重排序完成信号
    input  wire               inv,          // 0 表示FFT运算，1 表示IFFT运算
    input  wire signed [15:0] xm1_re,       // 输入数据 xm1 实部
    input  wire signed [15:0] xm1_im,       // 输入数据 xm1 虚部
    input  wire signed [15:0] xm2_re,       // 输入数据 xm2 实部
    input  wire signed [15:0] xm2_im,       // 输入数据 xm2 虚部
    input  wire signed [15:0] w_re,         // 旋转因子实部
    input  wire signed [15:0] w_im,         // 旋转因子虚部
    input  wire         [9:0] N,            // 点数
    input  wire         [3:0] M,            // 级数

// ============== 输出信号 ===============
    output wire         [8:0] rd_addr1,     // 输入数据 xm1 地址
    output wire         [8:0] rd_addr2,     // 输入数据 xm2 地址
    output wire         [8:0] wr_addr1,     // 输出数据 xn1 地址
    output wire         [8:0] wr_addr2,     // 输出数据 xn2 地址
    output wire         [7:0] rd_addr_wn,   // 旋转因子地址
    output wire               rd_en,        // 输入数据读取使能
    output wire               wr_en,        // 输出数据读取使能
    output wire               wn_rd_en,     // 旋转因子读取使能
    output wire signed [15:0] xn1_re,       // 输出数据 xn1 实部
    output wire signed [15:0] xn1_im,       // 输出数据 xn1 虚部
    output wire signed [15:0] xn2_re,       // 输出数据 xn2 实部
    output wire signed [15:0] xn2_im,       // 输出数据 xn2 虚部
    output reg                fft_start,
    output reg                fft_complete  // fft计算完成信号
    );

// ======================== 读写地址信号 ========================
    reg   [3:0] m;  // 级数计数器，最小为 1，最大为 9           1:M
    reg   [8:0] i;  // 组数计数器，最小为 1，最大为 256         1:2^(M-m)
    reg   [8:0] j;  // 蝶形在组内的编号，最小为 1，最大为 256    1:2^(m-1)
    wire [16:0] k;  // 上节点编号
    wire  [8:0] b;  // 下节点与上节点编号间隔
    wire [15:0] n;  // 旋转因子编号
    wire  [8:0] i_max;  // i 最大值
    wire  [8:0] j_max;  // j 最大值
    
    assign k = ((i-1) << m) + j;        // 未减 1
    assign b = (1 << (m-1));            // 未减 1
    assign n = (1 << (M-m)) * (j-1);    // 未减 1
    assign i_max = (1 << (M-m));
    assign j_max = (1 << (m-1));

    // 产生计算开始信号
    reg         cnt_mod2;   // 蝶形计算单元
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            fft_start <= 1'b0;
        else if (resort_complete == 1'b1)
            fft_start <= 1'b1;
        else if ((m == M) && (i == i_max) && (j == j_max) && (cnt_mod2 == 1'b1))
            fft_start <= 1'b0;
        else
            fft_start <= fft_start;
    end

    // 产生计算完成信号
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            fft_complete <= 1'b0;
        else if ((m == M) && (i == i_max) && (j == j_max) && (cnt_mod2 == 1'b1))
            fft_complete <= 1'b1;
        else
            fft_complete <= 1'b0;
    end

    // 模 2 计数器
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            cnt_mod2 <= 1'b0;
        else if (fft_start == 1'b1) begin
            if (cnt_mod2 == 1'b1)
                cnt_mod2 <= 1'b0;
            else
                cnt_mod2 <= cnt_mod2 + 1'b1;
        end else
            cnt_mod2 <= 1'b0;
    end

// ======================== 读写地址信号更新 ========================
    // m 信号更新
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            m <= 4'd1;
        else if ((fft_start == 1) && (i == i_max) && (j == j_max) && cnt_mod2 == 1'b1) begin
            if (m < M)
                m <= m + 1'b1;
            else
                m <= 4'd1;
        end else
            m <= m;
    end

    // i 信号更新
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            i <= 9'd1;
        else if ((fft_start == 1) && (j == j_max) && cnt_mod2 == 1'b1) begin
            if (i < i_max)
                i <= i + 1'b1;
            else
                i <= 9'd1;
        end else
            i <= i;
    end

    // j 信号更新
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            j <= 9'd1;
        else if ((fft_start == 1) && cnt_mod2 == 1'b1) begin
            if (j < j_max)
                j <= j + 1'b1;
            else
                j <= 9'd1;
        end else
            j <= j;
    end

    assign rd_en      = (fft_start == 1) ? (cnt_mod2 == 0 ? 1 : 0) : 1'b0;
    assign wr_en      = (fft_start == 1) ? (cnt_mod2 == 0 ? 0 : 1) : 1'b0;
    assign wn_rd_en   = (fft_start == 1) ? 1 : 1'b0;

    assign rd_addr1   = (fft_start == 1) ? k-1   : 9'dz;
    assign rd_addr2   = (fft_start == 1) ? k+b-1 : 9'dz;
    assign wr_addr1   = (fft_start == 1) ? k-1   : 9'dz;
    assign wr_addr2   = (fft_start == 1) ? k+b-1 : 9'dz;
    assign rd_addr_wn = (fft_start == 1) ? (n <<< (9 - M)) : 8'dz;

// ======================== 蝶形运算 ========================
    wire               butterfly_vld;

    butterfly u_butterfly (
        .clk    (clk),
        .rst_n  (rst_n),
        .en     (fft_start),
        .x1_re  (xm1_re),
        .x1_im  (xm1_im),
        .x2_re  (xm2_re),
        .x2_im  (xm2_im),
        .w_re   (w_re),
        .w_im   (w_im),
        .vld    (butterfly_vld),
        .y1_re  (xn1_re),
        .y1_im  (xn1_im),
        .y2_re  (xn2_re),
        .y2_im  (xn2_im) 
    );

endmodule
