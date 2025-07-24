// ---------------------------------------------------------------------------
// Multi‑Cycle Pipelined DPD MAC Unit (Streaming Coeffs, Complex, with CORDIC)
// Computes one complex output y[n] per (M+1)*K cycles.
// Each term: a[m,k] × x[n-m] × |x[n-m]|^(k-1)
// Uses pipelined CORDIC module to compute magnitude accurately.
// ---------------------------------------------------------------------------

module dpd_mac_array #(
  parameter DATA_WIDTH  = 16,   // Width of real/imag samples (Q2.14)
  parameter COEFF_WIDTH = 16,   // Width of real/imag coeffs (Q2.14)
  parameter FRAC_SZ     = 12,   // Number of fractional bits
  parameter M           = 2,   
    parameter K           = 3 ,    // Polynomial order // Memory depth
  parameter TOTAL_TERMS = (M+1)*K,
  parameter TERM_IDX_WIDTH = (TOTAL_TERMS <= 1) ? 1 : 
                              (TOTAL_TERMS <= 2) ? 1 : 
                              (TOTAL_TERMS <= 4) ? 2 : 
                              (TOTAL_TERMS <= 8) ? 3 : 4

)(
  input  logic                             clk,
  input  logic                             rst_n,
  input  logic                             valid_in,
  input  logic signed [DATA_WIDTH-1:0]     sample_window_re [0:M],
  input  logic signed [DATA_WIDTH-1:0]     sample_window_im [0:M],
  input  logic signed [COEFF_WIDTH-1:0]    coeff_in_re,
  input  logic signed [COEFF_WIDTH-1:0]    coeff_in_im,

  output logic                             coeff_req,
  output logic                             mac_busy,
  output logic                             valid_out,
  output logic signed [DATA_WIDTH-1:0]     y_dpd_re,
  output logic signed [DATA_WIDTH-1:0]     y_dpd_im,
  output logic [TERM_IDX_WIDTH-1:0]        term_idx  // Fixed: moved after width calculation
);

  // Calculate bit widths using parameters
 
  
  // Calculate bit widths with constant expressions
  localparam M_WIDTH = (M+1 <= 1) ? 1 : 
                       (M+1 <= 2) ? 1 : 
                       (M+1 <= 4) ? 2 : 
                       (M+1 <= 8) ? 3 : 4;
                       
  localparam K_WIDTH = (K+1 <= 1) ? 1 : 
                       (K+1 <= 2) ? 1 : 
                       (K+1 <= 4) ? 2 : 
                       (K+1 <= 8) ? 3 : 4;
                       


  // FSM states
  typedef enum logic [3:0] {
    IDLE,
    NEW_INPUT,
    WAIT_COEFF,
    LAUNCH_CORDIC,
    WAIT_CORDIC,
    POWER_INIT,
    POWER_NEXT,
    BASIS,
    MUL,
    ACC,
    FINISH
  } state_t;

  state_t state;

  // track whether we've computed magnitude for the current m
  logic m_calc_done;

  // Term-to-(m,k) LUTs with calculated widths
  logic [M_WIDTH-1:0] m_lut [0:TOTAL_TERMS-1];
  logic [K_WIDTH-1:0] k_lut [0:TOTAL_TERMS-1];

  // Initialize LUTs
  initial begin
    for (int t = 0; t < TOTAL_TERMS; t++) begin
      m_lut[t] = t / K;         // memory tap index
      k_lut[t] = (t % K) + 1;   // polynomial power
    end
  end

  // Request next coefficient in WAIT_COEFF
  assign coeff_req = (state == WAIT_COEFF);

  // Data-path regs with calculated widths
  logic [M_WIDTH-1:0]         m;
  logic [K_WIDTH-1:0]         k;
  logic signed [DATA_WIDTH-1:0]   x_m_re, x_m_im;
  logic signed [DATA_WIDTH+FRAC_SZ+1:0] magnitude;
  
  // Micro-pipeline registers
  logic [K_WIDTH-1:0]                     p_cnt;  // Fixed: K_WIDTH instead of $clog2(K)
  logic signed [DATA_WIDTH+FRAC_SZ+1:0]   power_reg;
  logic signed [DATA_WIDTH+FRAC_SZ+1:0]   basis_re_reg, basis_im_reg;
  logic signed [DATA_WIDTH+COEFF_WIDTH+FRAC_SZ+1:0] prod_re_reg, prod_im_reg;
  
  logic signed [DATA_WIDTH+COEFF_WIDTH+FRAC_SZ+1:0] acc_re, acc_im;

  // CORDIC interface
  logic                         cordic_valid_in, cordic_valid_out;
  logic signed [DATA_WIDTH-1:0] cordic_re_in, cordic_im_in;
  logic signed [DATA_WIDTH-1:0] cordic_mag_out;

  assign cordic_valid_in = (state == LAUNCH_CORDIC) && !cordic_valid_out;

  // latch inputs to CORDIC
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      cordic_re_in <= 0;
      cordic_im_in <= 0;
    end else if (state == WAIT_COEFF) begin
      cordic_re_in <= x_m_re;
      cordic_im_in <= x_m_im;
    end
  end

  cordic_mag #(
    .DATA_WIDTH(DATA_WIDTH),
    .ITER(16),
    .FRAC_SHIFT(FRAC_SZ)
  ) cordic_inst (
    .clk(clk),
    .rst_n(rst_n),
    .valid_in(cordic_valid_in),
    .x_in(cordic_re_in),
    .y_in(cordic_im_in),
    .valid_out(cordic_valid_out),
    .mag_out(cordic_mag_out)
  );

  // Saturation bounds
  logic signed [DATA_WIDTH-1:0] SAT_MAX = (1 <<< (DATA_WIDTH-1)) - 1;
  logic signed [DATA_WIDTH-1:0] SAT_MIN = -(1 <<< (DATA_WIDTH-1));

  // Main FSM + datapath
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      // reset everything
      state        <= IDLE;
      term_idx     <= 0;
      acc_re       <= 0;
      acc_im       <= 0;
      valid_out    <= 0;
      y_dpd_re     <= 0;
      y_dpd_im     <= 0;
      mac_busy     <= 0;
      magnitude    <= 0;
      x_m_re       <= 0;
      x_m_im       <= 0;
      m_calc_done  <= 0;
      p_cnt        <= 0;
      power_reg    <= 0;
      basis_re_reg <= 0;
      basis_im_reg <= 0;
      prod_re_reg  <= 0;
      prod_im_reg  <= 0;
    end else begin
      case (state)
        // -------------------------------------------------------------------
        IDLE: begin
          valid_out   <= 0;
          mac_busy    <= 0;
          if (valid_in) begin
            term_idx    <= 0;
            acc_re      <= 0;
            acc_im      <= 0;
            mac_busy    <= 1;
            m_calc_done <= 0;
            state       <= NEW_INPUT;
          end
        end

        NEW_INPUT: begin
          // latch tap order & sample
          m      <= m_lut[term_idx];
          k      <= k_lut[term_idx];
          x_m_re <= sample_window_re[m_lut[term_idx]];
          x_m_im <= sample_window_im[m_lut[term_idx]];
          state  <= WAIT_COEFF;
        end

        WAIT_COEFF: begin
          mac_busy <= 1;
          if (!m_calc_done)
            state <= LAUNCH_CORDIC;
          else
            state <= POWER_INIT;
        end

        LAUNCH_CORDIC: begin
          m_calc_done <= 1;
          state       <= WAIT_CORDIC;
        end

        WAIT_CORDIC: begin
          mac_busy <= 1;
          if (cordic_valid_out) begin
            magnitude <= cordic_mag_out;
            state     <= POWER_INIT;
          end
        end

        // -------------------------------------------------------------------
        // New micro‑pipeline to replace the old COMPUTE state
        POWER_INIT: begin
          // initialize power = |x|^0
          power_reg <= (1 <<< FRAC_SZ);
          p_cnt     <= 1;
          state     <= POWER_NEXT;
        end

        POWER_NEXT: begin
          // power_reg *= magnitude
          for (int p =1; p<k; p++)
          power_reg <= (power_reg * magnitude) >>> FRAC_SZ;
          // if (p_cnt == k-1) begin
          //   state <= BASIS;
          // end else begin
            p_cnt <= p_cnt + 1;
            state <= BASIS;
          //end
        end

        BASIS: begin
          // compute x * |x|^(k-1)
          basis_re_reg <= (x_m_re * power_reg) >>> FRAC_SZ;
          basis_im_reg <= (x_m_im * power_reg) >>> FRAC_SZ;
          state        <= MUL;
        end

        MUL: begin
          // multiply by coeff
          prod_re_reg <=  (basis_re_reg * coeff_in_re - basis_im_reg*coeff_in_im) >>> FRAC_SZ;
          prod_im_reg <= (basis_im_reg*coeff_in_re + basis_re_reg * coeff_in_im) >>> FRAC_SZ;
          state       <= ACC;
        end

        ACC: begin
          // accumulate
          acc_re <= acc_re + prod_re_reg;
          acc_im <= acc_im + prod_im_reg;
          if (k == K)
            m_calc_done <= 0;

          if (term_idx == TOTAL_TERMS-1) begin
            state <= FINISH;
          end else begin
            term_idx <= term_idx + 1;
            state    <= NEW_INPUT;
          end
        end

        // -------------------------------------------------------------------
        FINISH: begin
          mac_busy <= 0;
          // saturate and output
          y_dpd_re <= (acc_re > SAT_MAX ? SAT_MAX :
                       acc_re < SAT_MIN ? SAT_MIN : acc_re[DATA_WIDTH-1:0]);
          y_dpd_im <= (acc_im > SAT_MAX ? SAT_MAX :
                       acc_im < SAT_MIN ? SAT_MIN : acc_im[DATA_WIDTH-1:0]);
          valid_out <= 1;
          state     <= IDLE;
        end

      endcase
    end
  end

endmodule




// // -------------real numbers only-----------------------------
// // Multi-Cycle Pipelined DPD MAC Unit
// // Streams one coefficient per cycle
// // Uses LUTs for (m, k) mapping
// // ------------------------------------------

// module dpd_mac_array #(
//   parameter DATA_WIDTH  = 16,
//   parameter COEFF_WIDTH = 16,
//   parameter FRAC_SZ     = 14,
//   parameter M           = 1,
//   parameter K           = 3
// )(
//   input  logic                             clk,
//   input  logic                             rst_n,

//   input  logic                             valid_in,     // Start of new sample_window
//   input  logic signed [DATA_WIDTH-1:0]     sample_window [0:M],

//   input  logic                             coeff_valid,  // New coefficient available
//   input  logic signed [COEFF_WIDTH-1:0]    coeff_in,     // One coefficient per cycle

//   output logic                             coeff_req,    // Request next coefficient

//   output logic                             valid_out,
//   output logic signed [DATA_WIDTH-1:0]     y_dpd
// );

//   localparam TOTAL_TERMS = (M+1) * K;

//   // FSM state encoding
//   typedef enum logic [1:0] {
//     IDLE,
//     WAIT_COEFF,
//     COMPUTE,
//     FINISH
//   } state_t;

//   state_t state;

//   // Loop indexing
//   logic [$clog2(TOTAL_TERMS):0] term_idx;
//   logic [$clog2(M+1)-1:0] m;
//   logic [$clog2(K+1)-1:0] k;

//   // Lookup tables for m and k values
//   logic [$clog2(M+1)-1:0] m_lut [0:TOTAL_TERMS-1];
//   logic [$clog2(K+1)-1:0] k_lut [0:TOTAL_TERMS-1];

//   initial begin
//     for (int t = 0; t < TOTAL_TERMS; t++) begin
//       m_lut[t] = t / K;
//       k_lut[t] = (t % K) + 1;
//     end
//   end

//   // Internal registers
//   logic signed [DATA_WIDTH-1:0] x_m, abs_x;
//   logic signed [DATA_WIDTH+4:0] power;
//   logic signed [DATA_WIDTH+4:0] basis;
//   logic signed [DATA_WIDTH+COEFF_WIDTH+6:0] product;
//   logic signed [DATA_WIDTH+COEFF_WIDTH+8:0] acc;

//   assign coeff_req = (state == WAIT_COEFF);  // Pull one coeff per cycle

//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       state     <= IDLE;
//       term_idx  <= 0;
//       acc       <= 0;
//       y_dpd     <= 0;
//       valid_out <= 0;
//     end else begin
//       case (state)

//         IDLE: begin
//           valid_out <= 0;
//           if (valid_in) begin
//             term_idx <= 0;
//             acc      <= 0;
//             state    <= WAIT_COEFF;
//           end
//         end

//         WAIT_COEFF: begin
//           if (coeff_valid) begin
//             m = m_lut[term_idx];
//             k = k_lut[term_idx];

//             x_m   = sample_window[m];
//             abs_x = (x_m < 0) ? -x_m : x_m;

//             // |x|^(k-1)
//             power = 1;
//             for (int p = 0; p < (k - 1); p++)
//               power *= abs_x;

//             // basis = x * |x|^(k-1) >> FRAC_SZ
//             basis = (x_m * power) >>> FRAC_SZ;

//             // product = basis * coeff >> FRAC_SZ
//             product = (basis * coeff_in) >>> FRAC_SZ;

//             acc <= acc + product;

//             if (term_idx == TOTAL_TERMS - 1) begin
//               state <= FINISH;
//             end else begin
//               term_idx <= term_idx + 1;
//             end
//           end
//         end

//         FINISH: begin
//           // Saturate final output
//           logic signed [DATA_WIDTH-1:0] MAX = (1 <<< (DATA_WIDTH - 1)) - 1;
//           logic signed [DATA_WIDTH-1:0] MIN = -(1 <<< (DATA_WIDTH - 1));

//           if (acc > MAX) y_dpd <= MAX;
//           else if (acc < MIN) y_dpd <= MIN;
//           else y_dpd <= acc[DATA_WIDTH-1:0];

//           valid_out <= 1;
//           state     <= IDLE;
//         end

//         default: state <= IDLE;
//       endcase
//     end
//   end

// endmodule



// //////////////////////////buffering the coeffs first///////////////////////////////////
// module dpd_mac_array #(
//   parameter DATA_WIDTH  = 16,
//   parameter COEFF_WIDTH = 16,
//   parameter FRAC_SZ     = 14,
//   parameter M           = 1,
//   parameter K           = 3
// )(
//   input  logic                             clk,
//   input  logic                             rst_n,

//   input  logic                             valid_in,     // Input sample valid
//   input  logic signed [DATA_WIDTH-1:0]     sample_window [0:M],

//   input  logic                             coeff_valid,  // BRAM-style streaming
//   input  logic signed [COEFF_WIDTH-1:0]    coeff_in,     // One coefficient per cycle

//   output logic                             ready_for_coeffs, // High when accepting coeffs
//   output logic                             valid_out,
//   output logic signed [DATA_WIDTH-1:0]     y_dpd
// );

//   localparam TOTAL_TERMS = (M+1)*K;

//   // ----------------------------
//   // Internal state machine
//   // ----------------------------
//   typedef enum logic [1:0] {
//     IDLE,
//     LOAD_COEFFS,
//     COMPUTE,
//     FINISH
//   } state_t;

//   state_t state;

//   // ----------------------------
//   // Coefficient memory (register file)
//   // ----------------------------
//   logic signed [COEFF_WIDTH-1:0] coeffs_mem [0:TOTAL_TERMS-1];
//   logic [$clog2(TOTAL_TERMS):0]  coeff_idx;

//   // ----------------------------
//   // Term MAC computation
//   // ----------------------------
//   logic [$clog2(TOTAL_TERMS):0] term_idx;
//   logic [$clog2(M+1)-1:0]       m_lut [0:TOTAL_TERMS-1];
//   logic [$clog2(K+1)-1:0]       k_lut [0:TOTAL_TERMS-1];

//   initial begin
//     for (int t = 0; t < TOTAL_TERMS; t++) begin
//       m_lut[t] = t / K;
//       k_lut[t] = (t % K) + 1;
//     end
//   end

//   logic signed [DATA_WIDTH-1:0]           x_m, abs_x;
//   logic signed [DATA_WIDTH+4:0]           power;
//   logic signed [DATA_WIDTH+4:0]           basis;
//   logic signed [DATA_WIDTH+COEFF_WIDTH+6:0] product;
//   logic signed [DATA_WIDTH+COEFF_WIDTH+8:0] acc;

//   // ----------------------------
//   // FSM and datapath logic
//   // ----------------------------
//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       state          <= IDLE;
//       coeff_idx      <= 0;
//       term_idx       <= 0;
//       acc            <= 0;
//       y_dpd          <= 0;
//       valid_out      <= 0;
//       ready_for_coeffs <= 0;
//     end else begin
//       case (state)
//         // ------------------------
//         IDLE:
//         begin
//           valid_out <= 0;
//           if (valid_in) begin
//             coeff_idx <= 0;
//             ready_for_coeffs <= 1;
//             state <= LOAD_COEFFS;
//           end
//         end

//         // ------------------------
//         LOAD_COEFFS:
//         begin
//           if (coeff_valid) begin
//             coeffs_mem[coeff_idx] <= coeff_in;
//             coeff_idx <= coeff_idx + 1;

//             if (coeff_idx == TOTAL_TERMS - 1) begin
//               ready_for_coeffs <= 0;
//               term_idx <= 0;
//               acc <= 0;
//               state <= COMPUTE;
//             end
//           end
//         end

//         // ------------------------
//         COMPUTE:
//         begin
//           int k = k_lut[term_idx];
//           int m = m_lut[term_idx];

//           x_m   = sample_window[m];
//           abs_x = (x_m < 0) ? -x_m : x_m;

//           power = 1;
//           for (int p = 0; p < (k - 1); p++) power *= abs_x;

//           basis   = (x_m * power) >>> FRAC_SZ;
//           product = (basis * coeffs_mem[term_idx]) >>> FRAC_SZ;
//           acc    += product;

//           if (term_idx == TOTAL_TERMS - 1)
//             state <= FINISH;
//           else
//             term_idx <= term_idx + 1;
//         end

//         // ------------------------
//         FINISH:
//         begin
//           logic signed [DATA_WIDTH-1:0] MAX = (1 <<< (DATA_WIDTH - 1)) - 1;
//           logic signed [DATA_WIDTH-1:0] MIN = -(1 <<< (DATA_WIDTH - 1));

//           if (acc > MAX)
//             y_dpd <= MAX;
//           else if (acc < MIN)
//             y_dpd <= MIN;
//           else
//             y_dpd <= acc[DATA_WIDTH-1:0];

//           valid_out <= 1;
//           state     <= IDLE;
//         end

//         default: state <= IDLE;
//       endcase
//     end
//   end

// endmodule
