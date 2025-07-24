// ---------------------------------------------------------------------------
// Top module: single‑tap Sample Buffer → DPD MAC Array
//  • sample_buffer presents x[n−m] based on tap_idx
//  • dpd_mac_array consumes x_delayed_re/im, streams coeffs, produces y_dpd
// ---------------------------------------------------------------------------
module dpd_top #(
  parameter DATA_WIDTH  = 16,
  parameter COEFF_WIDTH = 16,
  parameter FRAC_SZ     = 12,
  parameter BUFF_WIDTH       = 4,
  parameter M           = 2,    // how many past taps we can select (0..M)
  parameter K           = 3     // polynomial order
)(
  input  logic                           clk,
  input  logic                           rst_n,
  input  logic                           ENABLE,

  // Incoming new sample x[n]
  input  logic                           buffer_valid_in,
  input  logic signed [DATA_WIDTH-1:0]   x_in_re,
  input  logic signed [DATA_WIDTH-1:0]   x_in_im,
  // // Coefficient streaming interface
  // input  logic                           coeff_valid,


  // Output interface
 
  output logic                           valid_out,
  output logic                            buf_full_re, buf_full_im, 
  //buf_empty_re, buf_empty_im,

  output logic signed [DATA_WIDTH-1:0]   y_dpd_re,
  output logic signed [DATA_WIDTH-1:0]   y_dpd_im
);
localparam TOTAL_TERMS =(M+1)*K;
logic mac_busy;
logic data_valid_re, data_valid_im;
logic [$clog2(TOTAL_TERMS)-1:0] term_idx;
logic signed [COEFF_WIDTH-1:0]  coeff_in_re;
logic signed [COEFF_WIDTH-1:0]  coeff_in_im;
logic coeff_req;

// ROM instance
 // logic signed [COEFF_WIDTH-1:0] coeff_r_re, coeff_r_im;
  coeff_rom #(
    .COEFF_WIDTH(COEFF_WIDTH),
    .DEPTH(TOTAL_TERMS)
  ) rom (
    .clk  (clk),
    .addr_re (term_idx),  
    .addr_im(term_idx),
    .enable (coeff_req),
    .data_re(coeff_in_re),
    .data_im(coeff_in_im)
  );

  // ----------------------------------------
  // Instantiate real‑path sample buffer
  // ----------------------------------------
  logic signed [DATA_WIDTH-1:0] x_delayed_re[0:M];
  sample_buffer #(
    .DATA_WIDTH(DATA_WIDTH),
    .WIDTH(BUFF_WIDTH)
  ) buf_re (
    .clk      (clk),
    .rst_n    (rst_n),
    .data_valid(data_valid_re),
    .mac_busy(mac_busy),
    .full(buf_full_re),
    //.empty(buf_empty_re),
    .valid_in (buffer_valid_in && ENABLE),
    .data_in  (x_in_re),
    .data_out (x_delayed_re)
  );

  // ----------------------------------------
  // Instantiate imag‑path sample buffer
  // ----------------------------------------
  logic signed [DATA_WIDTH-1:0] x_delayed_im[0:M];
  sample_buffer #(
    .DATA_WIDTH(DATA_WIDTH),
    .WIDTH(BUFF_WIDTH)
  ) buf_im (
    .clk      (clk),
    .rst_n    (rst_n),
    .full(buf_full_im),
   // .empty(buf_empty_im),
    .mac_busy(mac_busy),
    .data_valid(data_valid_im ),
    .valid_in (buffer_valid_in && ENABLE),
    .data_in  (x_in_im),
    .data_out (x_delayed_im)
  );
logic signed [DATA_WIDTH-1:0] latched_re[0:M];
logic signed [DATA_WIDTH-1:0] latched_im[0:M];
logic sample_valid_latched_re;
logic sample_valid_latched_im;


// always_ff @(posedge clk or negedge rst_n) begin
//   if (!rst_n) begin
//     sample_valid_latched_re <= 0;
//     sample_valid_latched_im <=0;
//       for (int i = 0; i < M; i++) begin
//       latched_re[i] <= '0;
//       latched_im[i] <='0;
//     end
//   end else begin
//     sample_valid_latched_re <= data_valid_re; // data_valid_re/im are in sync
//     sample_valid_latched_im <=data_valid_im;
//     latched_re <= x_delayed_re;
//     latched_im <= x_delayed_im;
//   end
// end
  // ----------------------------------------
  // Instantiate the simplified DPD MAC Array
  // ----------------------------------------
  dpd_mac_array #(
    .DATA_WIDTH (DATA_WIDTH),
    .COEFF_WIDTH(COEFF_WIDTH),
    .FRAC_SZ    (FRAC_SZ),
    .K          (K), 
    .M          (M)
  ) mac (
    .clk                (clk),
    .rst_n              (rst_n),
    .coeff_in_re        (coeff_in_re),
    .valid_in        ((data_valid_re || data_valid_im )&& ENABLE),
    .mac_busy        (mac_busy),
    .term_idx        (term_idx),
    .sample_window_re(x_delayed_re),
    .sample_window_im(x_delayed_im),
    .coeff_in_im        (coeff_in_im),
    .coeff_req          (coeff_req),     // pulses when MAC needs next coeff
    .valid_out          (valid_out),
    .y_dpd_re           (y_dpd_re),
    .y_dpd_im           (y_dpd_im)
  );

endmodule

