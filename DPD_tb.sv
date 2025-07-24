`timescale 1ns/1ps
module dpd_top_tb;

  // -------------------------------------------------------------------------
  // Parameters
  // -------------------------------------------------------------------------
  parameter DATA_WIDTH   = 16;
  parameter COEFF_WIDTH  = 16;
  parameter FRAC_SZ      = 12;
  parameter BUFF_WIDTH   = 3;       // buffer depth
  parameter M            = 2;
  parameter K            = 3;
  localparam TOTAL_TERMS = (M+1)*K;
  localparam CLK_PERIOD  = 10;      // 100 MHz
  localparam NUM_SAMPLES = 16;      // number of QAM points to stream

  // -------------------------------------------------------------------------
  // Signals
  // -------------------------------------------------------------------------
  logic                         clk, rst_n;
  logic                         ENABLE;
  logic                         buffer_valid_in;
  logic signed [DATA_WIDTH-1:0] x_in_re, x_in_im;
  logic                         valid_out;
  logic signed [DATA_WIDTH-1:0] y_dpd_re, y_dpd_im;
  logic mac_busy;
  logic                         buf_full_re, buf_full_im;

  // loop indices
  integer i;
  integer stream_idx;
  integer collect_idx;

  // sample & coefficient storage
  real sample_re  [0:NUM_SAMPLES-1];
  real sample_im  [0:NUM_SAMPLES-1];
  real coeffs_re_real [0:TOTAL_TERMS-1];
  real coeffs_im_real [0:TOTAL_TERMS-1];

  // File descriptors
  integer real_fd, imag_fd;

  // Memory images for fixed-point coeffs
  reg [COEFF_WIDTH-1:0] mem_re [0:TOTAL_TERMS-1];
  reg [COEFF_WIDTH-1:0] mem_im [0:TOTAL_TERMS-1];

  real dut_r, dut_i, golden_re, golden_im;

  // -------------------------------------------------------------------------
  // DUT Instantiation
  // -------------------------------------------------------------------------
  dpd_top #(
    .DATA_WIDTH (DATA_WIDTH),
    .COEFF_WIDTH(COEFF_WIDTH),
    .FRAC_SZ    (FRAC_SZ),
    .BUFF_WIDTH (BUFF_WIDTH),
    .M          (M),
    .K          (K)
  ) dut (
    .clk              (clk),
    .rst_n            (rst_n),
    .ENABLE           (ENABLE),
    .buffer_valid_in  (buffer_valid_in),
    .x_in_re          (x_in_re),
    .x_in_im          (x_in_im),
  //  .mac_busy       (mac_busy),
    .valid_out        (valid_out),
    .y_dpd_re         (y_dpd_re),
    .y_dpd_im         (y_dpd_im),
    .buf_full_re      (buf_full_re),
    .buf_full_im      (buf_full_im)
  );

  // -------------------------------------------------------------------------
  // Clock & Reset
  // -------------------------------------------------------------------------
  initial begin
    clk = 0;
    forever #(CLK_PERIOD/2) clk = ~clk;
  end

  initial begin
    rst_n = 0;
    # (CLK_PERIOD);
    rst_n = 1;
  end

  // -------------------------------------------------------------------------
  // Golden-model compute
  // -------------------------------------------------------------------------
  function real abs_r(input real v);
    abs_r = (v < 0) ? -v : v;
  endfunction

  task compute_golden(
    input integer n,
    input real      samp_re  [0:NUM_SAMPLES-1],
    input real      samp_im  [0:NUM_SAMPLES-1],
    input real      coeff_re [0:TOTAL_TERMS-1],
    input real      coeff_im [0:TOTAL_TERMS-1],
    output real     out_re,
    output real     out_im
  );
    integer m, k, idx;
    real    xre, xim, mag, power, bas_re, bas_im;
    begin
      out_re = 0; out_im = 0;
      for (m = 0; m <= M; m++) begin
        if ((n - m) < 0) begin
          xre = 0.0; xim = 0.0;
        end else begin
          xre = samp_re[n - m];
          xim = samp_im[n - m];
        end
        mag = $sqrt(xre*xre + xim*xim);
        for (k = 1; k <= K; k++) begin
          idx   = m*K + (k-1);
          power = (k == 1) ? 1.0 : mag**(k-1);
          bas_re = xre * power;
          bas_im = xim * power;
          out_re += bas_re * coeff_re[idx] - bas_im * coeff_im[idx];
          out_im += bas_im * coeff_re[idx] + bas_re * coeff_im[idx];
        end
      end
    end
  endtask

  // -------------------------------------------------------------------------
  // Test sequence: load coeffs, stream samples, collect outputs
  // -------------------------------------------------------------------------
  initial begin : TEST_SEQ
    // load reset and enable
    ENABLE = 0;
    buffer_valid_in = 0;
    # (2*CLK_PERIOD);

    // read and convert coefficients
    $display("Reading real coeffs...");
    $readmemh("coeffs_real.mem", mem_re);
    $display("Reading imag coeffs...");
    $readmemh("coeffs_imag.mem", mem_im);
    for (i = 0; i < TOTAL_TERMS; i++) begin
      coeffs_re_real[i] = $itor($signed(mem_re[i])) / (2.0**FRAC_SZ);
      coeffs_im_real[i] = $itor($signed(mem_im[i])) / (2.0**FRAC_SZ);
    end

    // // define QAM points
    // sample_re = '{ -3.0, -3.0, -3.0, -3.0,
    //                -1.0, -1.0, -1.0, -1.0,
    //                 1.0,  1.0,  1.0,  1.0,
    //                 3.0,  3.0,  3.0,  3.0 };
    // sample_im = '{ -3.0, -1.0,  1.0,  3.0,
    //                -3.0, -1.0,  1.0,  3.0,
    //                -3.0, -1.0,  1.0,  3.0,
    //                -3.0, -1.0,  1.0,  3.0 };

        // define QAM points
    sample_re = '{ -3.0, -3.0, -3.0, -3.0,
                   -1.0, -1.0, -1.0, -1.0,
                    1.0,  1.0,  1.0,  1.0,
                    3.0,  3.0,  3.0,  3.0 };
    sample_im = '{ -3.0, -1.0,  1.0,  3.0,
                   -3.0, -1.0,  1.0,  3.0,
                   -3.0, -1.0,  1.0,  3.0,
                   -3.0, -1.0,  1.0,  3.0 };

    // open files
    real_fd = $fopen("dpd_real.txt", "w");
    imag_fd = $fopen("dpd_imag.txt", "w");
    if (real_fd == 0 || imag_fd == 0) $error("File open failed");

    // Phase 1: stream samples
    ENABLE = 1;
    buffer_valid_in = 1;
    stream_idx = 0;

     for (int L=0; L<BUFF_WIDTH;L++) begin
        x_in_re = $rtoi(sample_re[L] * (2.0**FRAC_SZ));
        x_in_im = $rtoi(sample_im[L] * (2.0**FRAC_SZ));
        @(posedge clk);
      end

    while (stream_idx < NUM_SAMPLES ) begin

      if(!buf_full_re && !buf_full_im ) begin
        // buffer_valid_in = 1;

        x_in_re = $rtoi(sample_re[stream_idx+BUFF_WIDTH] * (2.0**FRAC_SZ));
        x_in_im = $rtoi(sample_im[stream_idx+BUFF_WIDTH] * (2.0**FRAC_SZ));

        wait (valid_out);
        

        compute_golden(stream_idx,
                       sample_re, sample_im,
                       coeffs_re_real, coeffs_im_real,
                       golden_re, golden_im);

        dut_r = $itor(y_dpd_re) / (2.0**FRAC_SZ);
        dut_i = $itor(y_dpd_im) / (2.0**FRAC_SZ);

        $display("Sample %0d: DUT=(%0f,%0f) Golden=(%0f,%0f)",
                 stream_idx, dut_r, dut_i, golden_re, golden_im);
        $fdisplay(real_fd, "%0f", dut_r);
        $fdisplay(imag_fd, "%0f", dut_i);
        if (abs_r(golden_re - dut_r) > 0.5e-0 || abs_r(golden_im - dut_i) > 0.5e-0)
          $error("Mismatch at index %0d", stream_idx);
        else
          $display("PASS index %0d", stream_idx);
        stream_idx++;
      end
      @(posedge clk);
    end
    buffer_valid_in = 1;
    #10;
    buffer_valid_in = 0;

    // Phase 2: collect outputs
    collect_idx = 0;
     for (i = 0; i < NUM_SAMPLES; i++) begin
      @(posedge clk);
      

        collect_idx++;
        if (collect_idx == NUM_SAMPLES) begin
          $display("All %0d samples done.", NUM_SAMPLES);
          $fclose(real_fd);
          $fclose(imag_fd);
          $stop;
        end
    end
  end

endmodule

// works but only one sample
//`timescale 1ns/1ps
// module dpd_top_tb;

//   // -------------------------------------------------------------------------
//   // Parameters (keep in sync with dpd_top instantiation)
//   // -------------------------------------------------------------------------
//   parameter DATA_WIDTH   = 16;
//   parameter COEFF_WIDTH  = 16;
//   parameter FRAC_SZ      = 12;
//   parameter BUFF_WIDTH   = 4;
//   parameter M            = 1;
//   parameter K            = 3;
//   localparam TOTAL_TERMS = (M+1)*K;
//   localparam CLK_PERIOD  = 10;  // 100 MHz

//   // -------------------------------------------------------------------------
//   // Signals
//   // -------------------------------------------------------------------------
//   logic                           clk, rst_n;
//   logic                           ENABLE;
//   logic                           buffer_valid_in;
//   logic signed [DATA_WIDTH-1:0]   x_in_re, x_in_im;
//   logic                           coeff_valid;
//   logic signed [COEFF_WIDTH-1:0]  coeff_in_re, coeff_in_im;
//   logic                           coeff_req;
//   logic                           valid_out;
//   logic                           mac_busy;
//   logic signed [DATA_WIDTH-1:0]   y_dpd_re, y_dpd_im;
//   real dut_re;
//   real dut_im;

//   // -------------------------------------------------------------------------
//   // DUT Instantiation
//   // -------------------------------------------------------------------------
//   dpd_top #(
//     .DATA_WIDTH (DATA_WIDTH),
//     .COEFF_WIDTH(COEFF_WIDTH),
//     .FRAC_SZ    (FRAC_SZ),
//     .BUFF_WIDTH (BUFF_WIDTH),
//     .M          (M),
//     .K          (K)
//   ) dut (
//     .clk              (clk),
//     .rst_n            (rst_n),
//     .ENABLE           (ENABLE),
//     .buffer_valid_in  (buffer_valid_in),
//     .x_in_re          (x_in_re),
//     .x_in_im          (x_in_im),
//     .coeff_valid      (coeff_valid),
//     .coeff_in_re      (coeff_in_re),
//     .coeff_in_im      (coeff_in_im),
//     .coeff_req        (coeff_req),
//     .valid_out        (valid_out),
//     .mac_busy         (mac_busy),
//     .y_dpd_re         (y_dpd_re),
//     .y_dpd_im         (y_dpd_im)
//   );

//   // -------------------------------------------------------------------------
//   // Clock & Reset
//   // -------------------------------------------------------------------------
//   initial begin
//     clk = 0;
//     forever #(CLK_PERIOD/2) clk = ~clk;
//   end

//   initial begin
//     // Assert reset
//     rst_n = 0;
//     ENABLE = 0;
//     buffer_valid_in = 0;
//     coeff_valid = 0;
//     x_in_re = 0; x_in_im = 0;
//     coeff_in_re = 0; coeff_in_im = 0;
//     # (CLK_PERIOD * 5);

//     // Release reset and enable
//     rst_n = 1;
//     ENABLE = 1;
//   end

//   // -------------------------------------------------------------------------
//   // Golden‐model for single input case
//   // -------------------------------------------------------------------------
//   real golden_re, golden_im;
//   task compute_golden(
//     input real xre, xim,
//     input real coeffs_re [0:TOTAL_TERMS-1],
//     input real coeffs_im [0:TOTAL_TERMS-1],
//     output real out_re, out_im
//   );
//     real mag, power, bas_re, bas_im;
//     integer idx, m_idx, k_idx;
//     begin
//       out_re = 0;
//       out_im = 0;
//       // Only one valid sample inserted at tap index 0; older taps are zero
//       for (idx = 0; idx < TOTAL_TERMS; idx++) begin
//         m_idx = idx / K;       // tap index
//         k_idx = (idx % K) + 1; // power
//         // only nonzero when m_idx == 0
//         if (m_idx == 0) begin
//           mag     = $sqrt(xre*xre + xim*xim);
//           power   = (k_idx == 1) ? mag :
//                     (k_idx == 2) ? mag*mag :
//                                    mag*mag*mag;
//           bas_re  = xre * power;
//           bas_im  = xim * power;
//           out_re += bas_re * coeffs_re[idx] - bas_im * coeffs_im[idx];
//           out_im += bas_re * coeffs_im[idx] + bas_im * coeffs_re[idx];
//         end
//       end
//     end
//   endtask

//   // -------------------------------------------------------------------------
//   // Main Test Sequence
//   // -------------------------------------------------------------------------
//   initial begin : TEST_SEQ
//     integer i;
//     real   coeffs_re_real [0:TOTAL_TERMS-1];
//     real   coeffs_im_real [0:TOTAL_TERMS-1];

//     // 1) Drive one complex sample: x = 1.0 + j*0.5
//     buffer_valid_in = 1;
//     x_in_re <= $rtoi(1.0 * (2.0**FRAC_SZ));
//     x_in_im <= $rtoi(0.5 * (2.0**FRAC_SZ));
//     @(posedge clk);
    

//     // 2) Provide a simple coefficient set when requested
//     //    e.g. coeff[i] = 1.0 for real, 0.0 for imag
//     for (i = 0; i < TOTAL_TERMS; i++) begin
//       coeffs_re_real[i] = 1.0;
//       coeffs_im_real[i] = 0.0;
//     end

//     // 3) Stream coefficients on each coeff_req
//     i = 0;
//     wait (coeff_req);
//     forever begin
//       if (coeff_req) begin
//         coeff_valid <= 1;
//         coeff_in_re <= $rtoi(coeffs_re_real[i] * (2.0**FRAC_SZ));
//         coeff_in_im <= $rtoi(coeffs_im_real[i] * (2.0**FRAC_SZ));
//         @(posedge clk);
//         coeff_valid <= 0;
//         i = i + 1;
//       end
//       if (i == TOTAL_TERMS) break;
//       @(posedge clk);
//     end
//   buffer_valid_in = 0;
//     // 4) Wait for valid_out
//     wait (valid_out);
//     @(posedge clk);

//     // 5) Compute golden output and compare
//     compute_golden(
//       1.0, 0.5,
//       coeffs_re_real,
//       coeffs_im_real,
//       golden_re,
//       golden_im
//     );

//     // convert DUT output back to real
//     dut_re = $itor(y_dpd_re) / (2.0**FRAC_SZ);
//     dut_im = $itor(y_dpd_im) / (2.0**FRAC_SZ);

//     $display("Golden = (%0f, %0f), DUT = (%0f, %0f)",
//              golden_re, golden_im, dut_re, dut_im);

//     if (abs(golden_re - dut_re) > 1.0/(2.0**FRAC_SZ) ||
//         abs(golden_im - dut_im) > 1.0/(2.0**FRAC_SZ)) begin
//       $error("Mismatch!");
//     end else begin
//       $display("PASS");
//     end

//     $stop;
//   end

//   // simple real abs function
//   function real abs(input real v);
//     abs = v < 0 ? -v : v;
//   endfunction

// endmodule





// // Enhanced self‑checking testbench for dpd_top with non‑zero fixed‑point scenario
// `timescale 1ns/1ps

// module dpd_top_tb;
//   //----------------------------------------------------------------------
//   // Parameters (must match dpd_top instantiation)
//   parameter DATA_WIDTH  = 16;
//   parameter COEFF_WIDTH = 16;
//   parameter FRAC_SZ     = 12;
//   parameter M           = 1;
//   parameter K           = 3;
//   localparam CLK_PERIOD = 10;

//   // Fixed‑point constants (Q2.14)
//   localparam logic signed [DATA_WIDTH-1:0] FP_ONE   = 1 << FRAC_SZ;      // 1.0
//   localparam logic signed [DATA_WIDTH-1:0] FP_HALF  = FP_ONE >>> 1;       // 0.5
//   localparam logic signed [DATA_WIDTH-1:0] FP_NEG1  = -FP_ONE;            // -1.0
//   // Expected output: y = (-1)*3 + (0.5 + 0.25 + 0.125) + (1 + 1 + 1) = 0.875
//   localparam logic signed [DATA_WIDTH-1:0] EXP_Y    = (FP_ONE * 7) / 8;    // 0.875 * 2^14 = 14336

//   //----------------------------------------------------------------------
//   // Testbench signals
//   logic clk;
//   logic rst_n;
//   logic ENABLE;
//   logic                 buffer_valid_in;
//   logic signed [DATA_WIDTH-1:0] x_in_re;
//   logic signed [DATA_WIDTH-1:0] x_in_im;

//   logic                 coeff_valid;
//   logic signed [COEFF_WIDTH-1:0] coeff_in_re;
//   logic signed [COEFF_WIDTH-1:0] coeff_in_im;

//   logic                 coeff_req;
//   logic                 valid_out;
//   logic signed [DATA_WIDTH-1:0] y_dpd_re;
//   logic signed [DATA_WIDTH-1:0] y_dpd_im;
//   int coeff_count;
// logic mac_busy;
//   //----------------------------------------------------------------------
//   // Custom abs function (not using $abs)
//   function automatic logic signed [DATA_WIDTH-1:0] abs_val(input logic signed [DATA_WIDTH-1:0] v);
//     return (v < 0) ? -v : v;
//   endfunction

//   //----------------------------------------------------------------------
//   // Instantiate DUT
//   dpd_top #(
//     .DATA_WIDTH(DATA_WIDTH),
//     .COEFF_WIDTH(COEFF_WIDTH),
//     .FRAC_SZ(FRAC_SZ),
//     .M(M),
//     .K(K)
//   ) dut (
//     .clk            (clk),
//     .rst_n          (rst_n),
//     .ENABLE         (ENABLE),
//     .buffer_valid_in(buffer_valid_in),
//     .x_in_re        (x_in_re),
//     .x_in_im        (x_in_im),
//     .coeff_valid    (coeff_valid),
//     .mac_busy       (mac_busy),
//     .coeff_in_re    (coeff_in_re),
//     .coeff_in_im    (coeff_in_im),
//     .coeff_req      (coeff_req),
//     .valid_out      (valid_out),
//     .y_dpd_re       (y_dpd_re),
//     .y_dpd_im       (y_dpd_im)
//   );

//   //----------------------------------------------------------------------
//   // Clock generation
//   initial begin
//     clk = 0;
//     forever #(CLK_PERIOD/2) clk = ~clk;
//   end

//   //----------------------------------------------------------------------
//   // Reset sequence
//   initial begin
//     rst_n           = 0;
//     buffer_valid_in = 0;
//     coeff_valid     = 0;
//     x_in_re         = '0;
//     x_in_im         = '0;
//     coeff_in_re     = '0;
//     coeff_in_im     = '0;
//     ENABLE          ='0;
//    @(posedge clk);
//     rst_n = 1;
//     ENABLE          ='1;
//   end

//   //----------------------------------------------------------------------
//   // Stimulus & self‑checking
//   initial begin
 
//     // 1) Feed three fixed‑point samples into sample_buffer
//     //    First: 1.0, then 0.5, then -1.0
//     @(posedge clk);
//     buffer_valid_in = 1; x_in_re = FP_ONE;  x_in_im = FP_ONE;
//      while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_ONE;
//         coeff_count++;
//         // #300;

//       end
//     end
  
//    // wait(!mac_busy)
//    // #300;
//     buffer_valid_in = 1; x_in_re = FP_HALF; x_in_im = 0;
//     //@(posedge clk);
//          while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_NEG1;
//         coeff_count++;
      
//       end
//     end

//     // wait(!mac_busy)
//     // #300;
//     buffer_valid_in = 1; x_in_re = FP_NEG1; x_in_im = FP_ONE;

//       while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_NEG1;
//         coeff_count++;
//         // #300;

//       end
//     end

//     // wait(!mac_busy)
//     // #300;
//     buffer_valid_in = 1; x_in_re = 3<<FRAC_SZ; x_in_im = 0;
//     //@(posedge clk);

//       while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_NEG1;
//         coeff_count++;
//         // #300;

//       end
//     end

//     // wait(!mac_busy)
//     // #300;
//     buffer_valid_in = 1; x_in_re =FP_NEG1; x_in_im = FP_ONE;
//     //@(posedge clk);
//     // @(posedge clk);
//     // buffer_valid_in = 1; x_in_re = FP_NEG1; x_in_im = 0;
//     //#300;    
//     //@(posedge clk);
//       while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_NEG1;
//         coeff_count++;
//          #300;

//       end
//     end

//   ENABLE          ='0;
//     // Deassert valid
//     @(posedge clk);
//     buffer_valid_in = 0; x_in_re = 0; x_in_im = 0;

//     // 2) Provide coefficients (all = 1.0) when requested
//     //    Total terms = (M+1)*K = 9
//     coeff_count = 0;
//     while (!valid_out) begin
//       @(posedge clk);
//       if (coeff_req) begin
//         coeff_valid = 1;
//         coeff_in_re = FP_ONE;
//         coeff_in_im = FP_NEG1;
//         coeff_count++;
//         #300;
//       end else begin
//         coeff_valid = 0;
//       end

//       if (coeff_count > ((M+1)*K + 2)*5) begin
//         $error("Stuck waiting for coefficients");
//         $stop;
//       end
//     end

//     // 3) Check output one cycle after valid_out
//     @(posedge clk);
//     if (y_dpd_re !== EXP_Y || y_dpd_im !== 0) begin
//       $error("DPD output mismatch: got re=%0d, im=%0d; expected re=%0d, im=0", y_dpd_re, y_dpd_im, EXP_Y);
//     end else begin
//       $display("PASS: DPD output re=%0d, im=%0d", y_dpd_re, y_dpd_im);
//     end

//     $stop;
//   end

// endmodule
