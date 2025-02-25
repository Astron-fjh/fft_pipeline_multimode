`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: UCAS
// Engineer: Astron-fjh
// 
// Create Date: 2025/02/13 21:25:27
// Design Name: Butterfly Operation Module
// Module Name: butterfly
// Project Name: fft_pipeline_multimode
// Target Devices: 
// Tool Versions: 
// Description: 
//   This module implements the butterfly operation for the Radix-2 FFT algorithm.
//   It takes two complex input signals (x1 and x2), along with a rotation factor (w), 
//   and calculates the corresponding FFT outputs (y1 and y2).
//   The inputs and outputs are in 16-bit signed fixed-point representation. 
//   The results of the butterfly operation are calculated using the formula:
//   y1 = x1 + w * x2
//   y2 = x1 - w * x2
//   The module includes a clock enable signal (en), a reset signal (rst_n), and 
//   generates an output valid signal (vld).
//
// Dependencies: 
//   None
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module butterfly(
// ============== 时钟信号 ===============
    input  wire               clk,
    input  wire               rst_n,
    input  wire               en,

// ============== 输入信号 ===============
    input  wire signed [15:0] x1_re,    // xm1
    input  wire signed [15:0] x1_im,
    input  wire signed [15:0] x2_re,    // xm2
    input  wire signed [15:0] x2_im,
    input  wire signed [15:0] w_re,     // Wnr
    input  wire signed [15:0] w_im,

// ============== 输出信号 ===============
    output wire               vld,
    output wire signed [15:0] y1_re,    // xn1
    output wire signed [15:0] y1_im,
    output wire signed [15:0] y2_re,    // xn2
    output wire signed [15:0] y2_im
    );

    reg en_r;

    always @(posedge clk or negedge rst_n) begin                                        
        if(!rst_n)                               
            en_r <= 1'b0;                                
        else                      
            en_r <= en;                                     
    end
    
    // ============== 蝶形运算 ===============
    // 输入数据延迟一个周期，与权重数据对齐
    reg signed [15:0] x1_re_r, x1_im_r;
    reg signed [15:0] x2_re_r, x2_im_r;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            x1_re_r <= 16'd0;
            x1_im_r <= 16'd0;
            x2_re_r <= 16'd0;
            x2_im_r <= 16'd0;
        end else if (en) begin
            x1_re_r <= x1_re;
            x1_im_r <= x1_im;
            x2_re_r <= x2_re;
            x2_im_r <= x2_im;
        end
    end

    reg signed [31:0] y1_re_r, y1_im_r;
    reg signed [31:0] y2_re_r, y2_im_r;

    always @(*) begin
        if (!rst_n) begin
            y1_re_r = 32'd0;
            y1_im_r = 32'd0;
            y2_re_r = 32'd0;
            y2_im_r = 32'd0;
        end else if (en) begin
            y1_re_r = x1_re_r + ((x2_re_r * w_re) >>> 13) - ((x2_im_r * w_im) >>> 13);   // 由于 w 被放大了 2^13 倍
            y1_im_r = x1_im_r + ((x2_re_r * w_im) >>> 13) + ((x2_im_r * w_re) >>> 13);   // 因此需右移 13 位
            y2_re_r = x1_re_r - ((x2_re_r * w_re) >>> 13) + ((x2_im_r * w_im) >>> 13);
            y2_im_r = x1_im_r - ((x2_re_r * w_im) >>> 13) - ((x2_im_r * w_re) >>> 13);
        end
    end

    assign vld = en_r;
    assign y1_re = y1_re_r[15:0];
    assign y1_im = y1_im_r[15:0];
    assign y2_re = y2_re_r[15:0];
    assign y2_im = y2_im_r[15:0];
    
endmodule
