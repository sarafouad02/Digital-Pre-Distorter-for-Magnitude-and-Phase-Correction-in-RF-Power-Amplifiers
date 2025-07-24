/*
  FSM-based sequential CORDIC magnitude calculator
  Maintains identical functionality to the pipelined version,
  but uses an iteration-counter FSM instead of unrolled pipeline stages.
*/
module cordic_mag #(
  parameter DATA_WIDTH = 16,   // input fixed-point width Q(format)
  parameter ITER       = 16,   // number of CORDIC iterations
  parameter FRAC_SHIFT = 12    // shift for gain normalization
)(
  input  logic                     clk,
  input  logic                     rst_n,
  input  logic                     valid_in,
  input  logic signed [DATA_WIDTH-1:0] x_in,
  input  logic signed [DATA_WIDTH-1:0] y_in,
  output logic                     valid_out,
  output logic signed [DATA_WIDTH-1:0] mag_out
);

  // internal widths
  localparam INT_W   = 5;
  localparam FULL_W  = 1 + INT_W + FRAC_SHIFT;
  localparam logic signed [15:0] RECIP_K = 16'sh09B7; // Q12

  // FSM states
  typedef enum logic [1:0] {
    STATE_IDLE,
    STATE_ROTATE,
    STATE_NORMALIZE
  } state_t;
  state_t state, next_state;

  // registers for iteration, data path
  logic [$clog2(ITER):0] iter_cnt;
  logic signed [FULL_W-1:0] x_reg, y_reg;
  logic signed [31:0] mult_reg;
  logic out_ready;

  // FSM sequential
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      state    <= STATE_IDLE;
      iter_cnt <= '0;
      x_reg    <= '0;
      y_reg    <= '0;
      mult_reg <= '0;
      out_ready<= 1'b0;
     // start <=1'b0;
    end else begin
      state    <= next_state;
      case (state)
        STATE_IDLE: begin
          out_ready <= 1'b0;
          //start <=1'b1;
          if (valid_in) begin
            // load absolute inputs
            x_reg <= (x_in < 0) ? -x_in : x_in;
            y_reg <= (y_in < 0) ? -y_in : y_in;
            iter_cnt <= '0;
          end
        end
        STATE_ROTATE: begin
        	//start<=0;
          // perform one CORDIC rotation
          if (y_reg >= 0) begin
            x_reg <= x_reg + (y_reg >>> iter_cnt);
            y_reg <= y_reg - (x_reg >>> iter_cnt);
          end else begin
            x_reg <= x_reg - (y_reg >>> iter_cnt);
            y_reg <= y_reg + (x_reg >>> iter_cnt);
          end
          iter_cnt <= iter_cnt + 1;
        end
        STATE_NORMALIZE: begin
          // apply gain correction
          mult_reg <= (x_reg * RECIP_K) >>> FRAC_SHIFT;
          out_ready<= 1'b1;
        end
      endcase
    end
  end

  // FSM combinational
  always_comb begin
    next_state = state;
    case (state)
      STATE_IDLE: if (valid_in)  next_state = STATE_ROTATE;
      STATE_ROTATE: if (iter_cnt == ITER) next_state = STATE_NORMALIZE;
      STATE_NORMALIZE: next_state = STATE_IDLE;
    endcase
  end

  // register output and clamp
  logic signed [DATA_WIDTH-1:0] mag_sat;
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      //mag_out   <= '0;
      valid_out <= 1'b0;
       mag_sat   <='0;
    end else begin
      valid_out <= out_ready;
      if (out_ready) begin
        // clamp to DATA_WIDTH signed
        if (mult_reg > ((1<<<(DATA_WIDTH-1))-1))
          mag_sat <= (1<<<(DATA_WIDTH-1))-1;
        else if (mult_reg < -((1<<<(DATA_WIDTH-1))))
          mag_sat <= -((1<<<(DATA_WIDTH-1)));
        else
          mag_sat <= mult_reg[DATA_WIDTH-1:0];
       // mag_out <= mag_sat;
      end
    end
  end
assign mag_out = mag_sat;
endmodule

// // ---------------------------------------------------------------------------
// // Pipelined CORDIC Magnitude Calculator Module
// // Computes |z| = sqrt(x^2 + y^2) via CORDIC vectoring mode
// // Latency = ITER + 1 cycles (pipelined), fixed-point arithmetic
// // ---------------------------------------------------------------------------
// module cordic_mag #(
//   parameter DATA_WIDTH = 16,   // Bit-width of input real/imag (signed fixed-point Q format)
//   parameter ITER       = 16,   // Number of CORDIC iterations (pipeline depth)
//   parameter FRAC_SHIFT = 14     // Right-shift to apply gain normalization
// )(
//   input  logic                           clk,       // Clock
//   input  logic                           rst_n,     // Active-low reset
//   input  logic                           valid_in,  // Input sample valid
//   input  logic signed [DATA_WIDTH-1:0]   x_in,      // Real part of input
//   input  logic signed [DATA_WIDTH-1:0]   y_in,      // Imag part of input
//   output logic                           valid_out, // Output valid after pipeline
//   output logic signed [DATA_WIDTH-1:0]   mag_out    // Computed magnitude |z|
// );

//   // ------------------------
//   // Pipeline registers
//   // ------------------------
//   localparam INT_W = 5;
//   localparam FULL_W = 1 + INT_W + FRAC_SHIFT;  // = 18

// logic signed [FULL_W-1:0] x_pipe [0:ITER];
// logic signed [FULL_W-1:0] y_pipe [0:ITER];

//   // logic signed [DATA_WIDTH:0] x_pipe   [0:ITER]; // the width needs to handle the number of iterations (shifts) to avoid overflow
//   // logic signed [DATA_WIDTH:0] y_pipe   [0:ITER];
//   logic                         valid_pipe[0:ITER];
//   logic signed [DATA_WIDTH-1:0] MAX = (1 <<< (DATA_WIDTH - 1)) - 1;  // Max positive value

//   // ------------------------
//   // Stage 0: Load inputs
//   // ------------------------
//   wire signed [FULL_W-1:0] x0 = (x_in < 0) ? -x_in : x_in;
//   wire signed [FULL_W-1:0] y0 = (y_in < 0) ? -y_in : y_in;

//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       x_pipe[0]     <= '0;
//       y_pipe[0]     <= '0;
//       valid_pipe[0] <= 1'b0;
//     end else begin
//       if (valid_in) begin
//         x_pipe[0]     <= x0;
//         y_pipe[0]     <= y0;
//         valid_pipe[0] <= 1'b1;
//       end else begin
//         valid_pipe[0] <= 1'b0;
//       end
//     end
//   end

//   // ------------------------
//   // CORDIC Iterations
//   // ------------------------

// generate
//   for (genvar i = 0; i < ITER; i++) begin : cordic_stages
//     always_ff @(posedge clk or negedge rst_n) begin
//       if (!rst_n) begin
//         x_pipe[i+1]     <= '0;
//         y_pipe[i+1]     <= '0;
//         valid_pipe[i+1] <= 1'b0;
//       end else begin
//         valid_pipe[i+1] <= valid_pipe[i];
//         if (valid_pipe[i]) begin
//           if (y_pipe[i] >= 0) begin
//             x_pipe[i+1] <= x_pipe[i] + (y_pipe[i] >>> i);
//             y_pipe[i+1] <= y_pipe[i] - (x_pipe[i] >>> i);
//           end else begin
//             x_pipe[i+1] <= x_pipe[i] - (y_pipe[i] >>> i);
//             y_pipe[i+1] <= y_pipe[i] + (x_pipe[i] >>> i);
//           end
//         end
//       end
//     end
//   end
// endgenerate

//   // ------------------------
//   // Final stage: Output magnitude
//   // ------------------------
//   localparam logic signed [15:0] RECIP_K = 16'sh09B7;  // Q12:
//    //16'sh26DD; //Q14
//   //16'sh09B7;  // Q12: 1/1.64676 ≈ 0.60725
//   logic signed [31:0] mult_pipe;
//   logic signed [31:0] norm_tmp;



// // Stage‑0 load x0, y0 instead of x_in, y_in.

//  logic                     valid_out_pipe;
//   // Multiply and delay valid
//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       mult_pipe       <= '0;
//       valid_out_pipe  <= 1'b0;
//     end else begin
//       valid_out_pipe <= valid_pipe[ITER];
//       if (valid_pipe[ITER]) begin
//         mult_pipe <= (x_pipe[ITER] * RECIP_K) >>> FRAC_SHIFT;
//       end
//     end
//   end

//   // Final register: clamp and assert valid_out only when data valid
//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       mag_out   <= '0;
//      valid_out <= 1'b0;
//     end else begin
//       if (valid_out_pipe) begin
//         // saturate
//         if (mult_pipe > MAX)
//           mag_out <= MAX;
//         else if (mult_pipe < -MAX)
//           mag_out <= -MAX;
//         else
//           mag_out <= mult_pipe[DATA_WIDTH-1:0];
       
//       end else begin
//        valid_out <= 1'b0;
//       end
//     end
//   end
// //assign valid_out = valid_out_pipe;

// endmodule