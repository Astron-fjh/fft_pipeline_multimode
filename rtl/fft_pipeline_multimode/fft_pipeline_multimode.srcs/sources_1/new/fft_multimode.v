`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: UCAS
// Engineer: Astron-fjh
// 
// Create Date: 2025/02/13 23:34:23
// Design Name: FFT Multimode Controller
// Module Name: fft_multimode
// Project Name: fft_pipeline_multimode
// Target Devices: 
// Tool Versions: 
// Description: 
//   This module implements a flexible FFT/IFFT (Fast Fourier Transform / Inverse Fast Fourier Transform) processor
//   supporting multiple modes with different point sizes. The module takes in complex input data and provides
//   the corresponding FFT or IFFT results. It includes data reordering, control signal generation, and 
//   communication with external memory for both input and output data storage. The module supports 
//   64-point, 128-point, 256-point, and 512-point FFT/IFFT operations based on the input control signals.
// 
// Dependencies: 
//   - SRAM for data storage
//   - Block RAM for storing rotation factors
//
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module fft_multimode(
// ============== 时钟信号 ===============
    input  wire               clk,
    input  wire               rst_n,

// ============== 输入信号 ===============
    input  wire               inv,         // 0 表示FFT运算，1 表示IFFT运算
    input  wire  [1:0]        np,          // FFT/IFFT点数：0表示64点；1表示128点；2表示256点；3表示512点
    input  wire               stb,         // 输入数据有效指示，低电平有效
    input  wire               sop_in,      // 每组输入数据第一个有效数据指示，高电平有效
    input  wire signed [15:0] x_re,        // 输入数据实部，二进制补码定点格式
    input  wire signed [15:0] x_im,        // 输入数据虚部，二进制补码定点格式

// ============== 输出信号 ===============
    output wire               valid_out,   // 输出数据有效指示，低电平有效
    output wire               sop_out,     // 每组输出数据第一个有效数据指示，低电平有效
    output reg  signed [15:0] y_re,        // 输出数据实部，二进制补码定点格式
    output reg  signed [15:0] y_im         // 输出数据虚部，二进制补码定点格式
    );

    reg   [9:0] N;  // 点数
    reg   [3:0] M;  // 级数

    always @(*) begin
        case(np)
            2'd0: begin
                N = 10'd64;
                M = 4'd6;
            end
            2'd1: begin
                N = 10'd128;
                M = 4'd7;
            end
            2'd2: begin
                N = 10'd256;
                M = 4'd8;
            end
            2'd3: begin
                N = 10'd512;
                M = 4'd9;
            end
        endcase
    end

// ======================== 输入数据存储sram ========================
    wire               wr_en;
    wire               rd_en;
    wire         [8:0] wr_addr1;
    wire         [8:0] wr_addr2;
    wire         [8:0] rd_addr1;
    wire         [8:0] rd_addr2;
    wire signed [15:0] xn1_re_in;
    wire signed [15:0] xn1_im_in;
    wire signed [15:0] xn2_re_in;
    wire signed [15:0] xn2_im_in;

    wire signed [15:0] xm1_re_out;
    wire signed [15:0] xm1_im_out;
    wire signed [15:0] xm2_re_out;
    wire signed [15:0] xm2_im_out;

    sram_2x16x512bit sram_2x16x512bit_init(
        .clk        (clk),
        .rst_n      (rst_n),
        .wr_en      (wr_en),
        .rd_en      (rd_en),
        .wr_addr1   (wr_addr1),
        .wr_addr2   (wr_addr2),
        .rd_addr1   (rd_addr1),
        .rd_addr2   (rd_addr2),
        .xn1_re_in  (xn1_re_in),
        .xn1_im_in  (xn1_im_in),
        .xn2_re_in  (xn2_re_in),
        .xn2_im_in  (xn2_im_in),
        .xm1_re_out (xm1_re_out),
        .xm1_im_out (xm1_im_out),
        .xm2_re_out (xm2_re_out),
        .xm2_im_out (xm2_im_out)
    );

// ======================= 旋转因子存储 RAM =======================
    wire         [7:0] rd_addr_wn;
    wire signed [15:0] w_re, w_im;
    blk_mem_gen_0 u_w_re_ram (  // 存储 W_re 定点化数据
        .clka   (clk),          // input  wire          clka
        .wea    (1'b0),         // input  wire  [0 : 0] wea   高为写，低为读
        .addra  (),             // input  wire  [7 : 0] addra
        .dina   (),             // input  wire [15 : 0] dina
        .clkb   (clk),          // input  wire          clkb
        .addrb  (rd_addr_wn),   // input  wire  [7 : 0] addrb
        .doutb  (w_re)          // output wire [15 : 0] doutb
    );

    blk_mem_gen_1 u_w_im_ram (  // 存储 W_im 定点化数据
        .clka   (clk),          // input  wire          clka
        .wea    (1'b0),         // input  wire  [0 : 0] wea   高为写，低为读
        .addra  (),             // input  wire  [7 : 0] addra
        .dina   (),             // input  wire [15 : 0] dina
        .clkb   (clk),          // input  wire          clkb
        .addrb  (rd_addr_wn),   // input  wire  [7 : 0] addrb
        .doutb  (w_im)          // output wire [15 : 0] doutb
    );

// ========================== 数据重排 ==========================
    wire               resort_wr_en;
    wire               resort_complete;
    wire signed  [8:0] resort_wr_addr1;
    wire signed [15:0] resort_xn1_re_in, resort_xn1_im_in;
    in_resort u_in_resort (
        .clk                (clk),
        .rst_n              (rst_n),
        .N                  (N),
        .inv                (inv),
        .sop_in             (sop_in),
        .x_re               (x_re),
        .x_im               (x_im),
        .wr_en              (resort_wr_en),
        .resort_complete    (resort_complete),
        .addr_resort        (resort_wr_addr1),
        .x_re_resorted      (resort_xn1_re_in),
        .x_im_resorted      (resort_xn1_im_in)
    );

// ======================== FFT/IFFT 计算 ========================
    wire               fft_ctrl_rd_en,     fft_ctrl_wr_en;
    wire               fft_start,          fft_complete;
    wire         [8:0] fft_ctrl_rd_addr1,  fft_ctrl_wr_addr1;
    wire signed [15:0] fft_ctrl_xn1_re_in, fft_ctrl_xn1_im_in;
    fft_ctrl u_fft_ctrl(
        .clk                (clk),
        .rst_n              (rst_n),
        .resort_complete    (resort_complete),
        .inv                (inv),
        .xm1_re             (xm1_re_out),
        .xm1_im             (xm1_im_out),
        .xm2_re             (xm2_re_out),
        .xm2_im             (xm2_im_out),
        .w_re               (w_re),
        .w_im               (w_im),
        .N                  (N),
        .M                  (M),
        .rd_addr1           (fft_ctrl_rd_addr1),
        .rd_addr2           (rd_addr2),
        .wr_addr1           (fft_ctrl_wr_addr1),
        .wr_addr2           (wr_addr2),
        .rd_addr_wn         (rd_addr_wn),
        .rd_en              (fft_ctrl_rd_en),
        .wr_en              (fft_ctrl_wr_en),
        .wn_rd_en           (wn_rd_en),
        .xn1_re             (fft_ctrl_xn1_re_in),
        .xn1_im             (fft_ctrl_xn1_im_in),
        .xn2_re             (xn2_re_in),
        .xn2_im             (xn2_im_in),
        .fft_start          (fft_start),
        .fft_complete       (fft_complete)
    );

// ======================== 计算结果输出 ========================
    reg        dataout_start;
    reg  [9:0] cnt_dataout;
    wire       dataout_rd_en;
    wire [8:0] dataout_rd_addr1;

    always @(posedge clk or negedge rst_n) begin
        if (!rst_n)
            dataout_start <= 1'b0;
        else if (fft_complete == 1)
            dataout_start <= 1'b1;
        else if (cnt_dataout == N-1)
            dataout_start <= 1'b0;
        else
            dataout_start <= dataout_start;
    end

    always @(posedge clk or negedge rst_n) begin                                        
        if (!rst_n)                               
            cnt_dataout <= 10'd0;
        else if (dataout_start == 1'b1) begin
            if (cnt_dataout < N-1)
                cnt_dataout <= cnt_dataout + 1'b1;
            else
                cnt_dataout <= 10'd0;
        end else
            cnt_dataout <= 10'd0;           
    end

    assign dataout_rd_en    = dataout_start ? 1'b1 : 1'b0;
    assign dataout_rd_addr1 = dataout_start ? cnt_dataout : 9'dz;
    assign sop_out          = (dataout_start == 1 && cnt_dataout == 0) ? 1 : 0;

    assign valid_out = ~dataout_start;

    always @(*) begin
        y_re = 16'd0;
        y_im = 16'd0;
        if (dataout_start == 1'b1) begin
            if (inv == 0) begin
                y_re = xm1_re_out;
                y_im = xm1_im_out;
            end else begin
                case(np)
                    2'd0: begin
                        y_re = (xm1_re_out >>> 6);
                        y_im = ((-xm1_im_out) >>> 6);
                    end
                    2'd1: begin
                        y_re = (xm1_re_out >>> 7);
                        y_im = ((-xm1_im_out) >>> 7);
                    end
                    2'd2: begin
                        y_re = (xm1_re_out >>> 8);
                        y_im = ((-xm1_im_out) >>> 8);
                    end
                    2'd3: begin
                        y_re = (xm1_re_out >>> 9);
                        y_im = ((-xm1_im_out) >>> 9);
                    end
                endcase
            end
        end
    end

    assign rd_en    = dataout_start ? dataout_rd_en : fft_ctrl_rd_en;
    assign rd_addr1 = dataout_start ? dataout_rd_addr1 : fft_ctrl_rd_addr1;

    assign wr_en     = fft_start ? fft_ctrl_wr_en : resort_wr_en;
    assign wr_addr1  = fft_start ? fft_ctrl_wr_addr1 : resort_wr_addr1;
    assign xn1_re_in = fft_start ? fft_ctrl_xn1_re_in : resort_xn1_re_in;
    assign xn1_im_in = fft_start ? fft_ctrl_xn1_im_in : resort_xn1_im_in;

endmodule
