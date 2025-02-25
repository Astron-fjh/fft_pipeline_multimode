`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: UCAS
// Engineer: Astron-fjh
// 
// Create Date: 2025/02/14 23:38:59
// Design Name: SRAM Module for FFT Butterfly Operation
// Module Name: sram_2x16x512bit
// Project Name: fft_pipeline_multimode
// Target Devices: 
// Tool Versions: 
// Description: 
//   This module implements a dual-port SRAM with 512 memory locations, where each location 
//   stores 16-bit signed real and imaginary data separately. The SRAM is used to store the 
//   intermediate results of butterfly operations in FFT/IFFT processing. It supports both 
//   read and write operations with two independent address ports (wr_addr1, wr_addr2 for write, 
//   rd_addr1, rd_addr2 for read). The write operation occurs immediately with a one-cycle delay 
//   for each input data, while read operations output the stored data from the specified addresses.
//
//   The SRAM stores both real and imaginary parts of the butterfly operation outputs and 
//   provides them as inputs to the next stage of the FFT process. The data is read when 
//   `rd_en` is asserted and written when `wr_en` is asserted. 
//
// Dependencies: 
//   None
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module sram_2x16x512bit(
// ============== 时钟信号 ===============
    input  wire               clk,
    input  wire               rst_n,

// ============== 输入信号 ===============
    input  wire               wr_en,
    input  wire               rd_en,
    input  wire         [8:0] wr_addr1,
    input  wire         [8:0] wr_addr2,
    input  wire         [8:0] rd_addr1,
    input  wire         [8:0] rd_addr2,
    input  wire signed [15:0] xn1_re_in,   // 蝶形运算输出
    input  wire signed [15:0] xn1_im_in,
    input  wire signed [15:0] xn2_re_in,
    input  wire signed [15:0] xn2_im_in,

// ============== 输出信号 ===============
    output reg  signed [15:0] xm1_re_out,  // 蝶形运算输入
    output reg  signed [15:0] xm1_im_out,
    output reg  signed [15:0] xm2_re_out,
    output reg  signed [15:0] xm2_im_out
    );

    reg signed [15:0] re_sram[511:0];
    reg signed [15:0] im_sram[511:0];

    // write, delay for 1 cycle
    integer i;
    always@(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            for (i = 0; i < 512; i = i+1) begin
                re_sram[i] = 16'd0;
                im_sram[i] = 16'd0;
            end
        end else if (wr_en == 1) begin
            re_sram[wr_addr1] = xn1_re_in;
            im_sram[wr_addr1] = xn1_im_in;
            re_sram[wr_addr2] = xn2_re_in;
            im_sram[wr_addr2] = xn2_im_in;
        end
    end

    // read
    always @(*) begin
        if (rd_en == 1) begin
            xm1_re_out = re_sram[rd_addr1];
            xm1_im_out = im_sram[rd_addr1];
            xm2_re_out = re_sram[rd_addr2];
            xm2_im_out = im_sram[rd_addr2];
        end
    end
endmodule
