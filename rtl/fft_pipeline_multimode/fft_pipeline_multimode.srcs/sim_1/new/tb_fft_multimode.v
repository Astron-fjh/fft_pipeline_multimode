`timescale 1ns / 1ps
//////////////////////////////////////////////////////////////////////////////////
// Company: UCAS
// Engineer: Astron-fjh
// 
// Create Date: 2025/02/16 18:50:41
// Design Name: 
// Module Name: tb_fft_multimode
// Project Name: fft_pipeline_multimode
// Target Devices: 
// Tool Versions: 
// Description: 
//   This is the testbench for the fft_multimode module, which tests the FFT/IFFT
//   functionality with different point sizes (64, 128, 256, and 512 points). The 
//   testbench generates the necessary input signals and applies them to the DUT 
//   (Device Under Test). It also captures the output data (real and imaginary parts)
//   and writes them to external files for verification.
//
//   The testbench initializes the input data, applies the necessary control signals,
//   and verifies the output results for correctness. It also handles the reading and
//   writing of data to and from memory files, allowing for easy verification of 
//   the results in a simulation environment.
//
// Dependencies: 
//   - FFT/IFFT module (fft_multimode)
//   - Memory files (input and output data)
// 
// Revision:
// Revision 0.01 - File Created
// Additional Comments:
// 
//////////////////////////////////////////////////////////////////////////////////


module tb_fft_multimode;
    reg         clk;
    reg         rst_n;
    reg         inv;
    reg   [1:0] np;
    reg         stb;
    reg         sop_in;
    reg  [15:0] x_re;
    reg  [15:0] x_im;
    wire        valid_out;
    wire        sop_out;
    wire [15:0] y_re;
    wire [15:0] y_im;

    parameter  period = 200;

    reg  [9:0] N;
    reg [15:0] mem_re[511:0];
    reg [15:0] mem_im[511:0];

    fft_multimode u_fft_multimode(
        .clk        (clk),
        .rst_n      (rst_n),
        .inv        (inv),
        .np         (np),
        .stb        (stb),
        .sop_in     (sop_in),
        .x_re       (x_re),
        .x_im       (x_im),
        .valid_out  (valid_out),
        .sop_out    (sop_out),
        .y_re       (y_re),
        .y_im       (y_im)
    );

    integer fdyre, fdyim;
    integer cnt0;
    integer stop_flag;
    
    // 读入数据
    initial begin
        clk     = 1'b1;
        rst_n   = 1'b0;
        inv     = 1;        // FFT/IFFT
        np      = 2'b11;    // 64/128/256/512
        sop_in  = 1'b0;
        stb     = 1'b1;

        x_re    = 16'd0;
        x_im    = 16'd0;
        
        cnt0    = 0;
        stop_flag = 0;

        case(np)
            2'b00: N = 10'd64;
            2'b01: N = 10'd128;
            2'b10: N = 10'd256;
            2'b11: N = 10'd512;
        endcase

        $readmemh("E:/Vivado/fft_pipeline_multimode/fft_pipeline_multimode.srcs/sim_1/new/x_re512.txt", mem_re);
        $readmemh("E:/Vivado/fft_pipeline_multimode/fft_pipeline_multimode.srcs/sim_1/new/x_im512.txt", mem_im);

        #(period);
        rst_n = 1'b1;

        #(period);
        stb = 1'b0; // 低电平有效

        while(cnt0 < N) begin
            if(cnt0 == 0)
                sop_in = 1;
            else
                sop_in = 0;

            x_re = mem_re[cnt0];
            x_im = mem_im[cnt0];
            #(period);
            cnt0 = cnt0 + 1;
        end

        stb = 1'b1;
        stop_flag = 1;
    end

    // 写入计算结果
    integer cnt1;
    initial begin
        fdyre = $fopen("E:/Vivado/fft_pipeline_multimode/fft_pipeline_multimode.srcs/sim_1/new/y_ifft512_re.txt", "wb");
        fdyim = $fopen("E:/Vivado/fft_pipeline_multimode/fft_pipeline_multimode.srcs/sim_1/new/y_ifft512_im.txt", "wb");
        cnt1 = 0;
        
        while(1) begin
            if(valid_out == 0) begin
                $fwrite(fdyre, "%04x\n", y_re);
                $fwrite(fdyim, "%04x\n", y_im);
                cnt1 = cnt1 + 1;
            end else if(valid_out == 1 && stop_flag == 1 && cnt1 == N) begin
                $fclose(fdyre);
                $fclose(fdyim);
                $stop;
            end

            #(period);
        end
    end

    always #(period/2) clk = ~clk;

endmodule
