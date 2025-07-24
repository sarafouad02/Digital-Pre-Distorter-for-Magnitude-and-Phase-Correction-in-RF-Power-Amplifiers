// // ---------------------------------------------------------------------------
// // Parametric Sample Buffer with Tap Selection
// // Stores M+1 past samples and presents only the delayed sample indexed by tap_idx
// // ---------------------------------------------------------------------------
// module sample_buffer #(
//     parameter DATA_WIDTH = 16,  // bit‑width of the input samples
//     parameter M          = 1    // memory depth (number of past samples to store)
// )(
//     input  logic                             clk,       // clock
//     input  logic                             rst,     // async active‑low reset
//     input  logic                             valid_in,  // high when data_in is valid
//     input  logic signed [DATA_WIDTH-1:0]     data_in,   // new sample x[n]

//     input  logic [$clog2(M+1)-1:0]           tap_idx,   // which delayed tap to output (0..M)
//     output logic signed [DATA_WIDTH-1:0]     data_out   // selected delayed sample x[n-tap_idx]
// );

//   // Internal shift‑register storing x[n], x[n-1], …, x[n-M]
//   logic signed [DATA_WIDTH-1:0] buffer [M];

//   // -------------------------------------------------------------------------
//   // Shift logic: on valid_in, push in new sample at buffer[0], shift old samples
//   // -------------------------------------------------------------------------
//   always_ff @(posedge clk or negedge rst) begin
//         if (!rst) begin
//             for (int i = 0; i <= M; i++) begin
//                 buffer[i] <= '0;
//             end
//         end else if (valid_in) begin
//             // Shift older samples
//             for (int i = M; i > 0; i--) begin
//                 buffer[i] <= buffer[i-1];
//             end
//             // Store new input
//             buffer[0] <= data_in;
//         end
//     end

//   // -------------------------------------------------------------------------
//   // Tap selection: combinationally output the requested delayed sample
//   // -------------------------------------------------------------------------
//   always_comb begin
//     // Bound check (optional): if tap_idx > M, return zero
//     if (tap_idx <= M)
//       data_out = buffer[tap_idx];
//     else
//       data_out = '0;
//   end

// endmodule


// module sample_buffer #(
//   parameter DATA_WIDTH = 16,  // width of each sample
//   parameter WIDTH      = 4    // number of columns per row
// )(
//   input  logic                       clk,
//   input  logic                       rst_n,      // active‑low
//   input  logic                       valid_in,   // each pulse shifts-in one new sample
//   input  logic                       mac_busy,   // back‑pressure for streaming
//   input  logic signed [DATA_WIDTH-1:0] data_in,
//   output logic                       data_valid, // high when data_out is valid
//   output logic signed [DATA_WIDTH-1:0] data_out [0:1]  // [0]=row0[col], [1]=row1[col]
// );

//   // Storage: two rows × WIDTH columns
//   logic signed [DATA_WIDTH-1:0] row0 [WIDTH-1:0];
//   logic signed [DATA_WIDTH-1:0] row1 [WIDTH-1:0];

//   // read pointer and state flags
//   logic [$clog2(WIDTH)-1:0] rd_ptr;
//   logic                     loading;   // still shifting in
//   logic                     streaming; // ready to stream out
//   int c;

//   //------------------------------------------------------------------------
//   // 1) LOAD & SHIFT‑IN PHASE
//   //------------------------------------------------------------------------
//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       // clear everything on reset
//       for (int c = 0; c < WIDTH; c++) begin
//         row0[c] <= '0;
//         row1[c] <= '0;
//       end
//       rd_ptr    <= 0;
//       loading   <= 1'b1;
//       streaming <= 1'b0;
//     end else begin
//       if (valid_in && loading) begin
//         // vertical shift: row0 → row1
//         for (int c = 0; c < WIDTH; c++)
//           row0[c] <= row1[c];
//         // horizontal shift into row0
//         for (int c = WIDTH-1; c > 0; c--)
//           row1[c] <= row1[c-1];
//         row1[0] <= data_in;

//         // advance pointer; once we've loaded WIDTH samples, switch modes
//         rd_ptr  <= rd_ptr + 1;
//         if (rd_ptr == WIDTH-1) begin
//           loading   <= 1'b0;
//           streaming <= 1'b1;
//           rd_ptr    <= 0;
//         end
//       end
//     end
//   end

//   //------------------------------------------------------------------------
//   // 2) STREAM‑OUT PHASE
//   //------------------------------------------------------------------------
//   always_ff @(posedge clk or negedge rst_n) begin
//     if (!rst_n) begin
//       data_valid <= 1'b0;
//       rd_ptr     <= 0;
//     end else if (streaming && !mac_busy) begin
//       // present one column per cycle
//       data_valid <= 1'b1;
//       rd_ptr     <= rd_ptr + 1;

//       // once we've streamed all WIDTH columns, go back to loading
//       if (rd_ptr == WIDTH-1) begin
//         streaming <= 1'b0;
//         loading   <= 1'b1;
//         rd_ptr    <= 0;
//       end
//     end else begin
//       data_valid <= 1'b0;
//     end
//   end

//   //------------------------------------------------------------------------
//   // 3) DRIVE data_out
//   //------------------------------------------------------------------------
//   always_comb begin
//     if (streaming && data_valid && !mac_busy) begin
//       data_out[0] = row0[rd_ptr];
//       data_out[1] = row1[rd_ptr];
//     end else begin
//       data_out[0] = '0;
//       data_out[1] = '0;
//     end
//   end

// sample_buffer.v
// sample_buffer.v
`timescale 1ns/1ps
module sample_buffer #(
  parameter int DATA_WIDTH = 16,
  parameter int WIDTH      = 3,    // depth of circular buffer
  parameter int M          = 2     // number of �previous� samples
)(
  input  logic                        clk,
  input  logic                        rst_n,
  input  logic                        valid_in,
  input  logic                        mac_busy,
  input  logic signed [DATA_WIDTH-1:0] data_in,
  output logic                        data_valid,
  output logic                        full,           // high when buffer is full
 // high once initial fill complete
  // window of M+1 samples
  output logic signed [DATA_WIDTH-1:0] data_out [0:M]
);

  // -------------------------------------------------------------------
  // Storage: one buffer for each of the M+1 sample positions
  // -------------------------------------------------------------------
  logic signed [DATA_WIDTH-1:0] row_buf [0:M][0:WIDTH-1];
  logic signed [DATA_WIDTH-1:0] prev_samples [0:M-1];
  logic                        init_fill_done;
  logic [$clog2(WIDTH)-1:0] wr_ptr, rd_ptr;
  logic [$clog2(WIDTH):0]   count;
  logic                     read_issued;
  int inc, dec;
  // full flag
  logic full_reg, empty;
  assign full = (count >= WIDTH && !dec);
  assign empty = count ==0 && !inc;

logic init_fill_done_reg;
  // initial fill flag: asserted when count first reaches WIDTH
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      init_fill_done <= 0;
      init_fill_done_reg<=0;
  end
    else if (count >= WIDTH) begin
      init_fill_done <= 1;
      init_fill_done_reg<=init_fill_done;
    end
  end

always_ff @(posedge clk or negedge rst_n) begin : proc_full_reg
  if(~rst_n) begin
    full_reg <= 0;
  end else begin
    full_reg <= full ;
  end
end
  //------------------------------------------------------------------------
  // WRITE: on each valid_in & not full
  //------------------------------------------------------------------------

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      wr_ptr       <= 0;
      for (int i = 0; i < M; i++)
        prev_samples[i] <= '0;
    end else if (valid_in && (count < WIDTH) && (!full && !init_fill_done_reg || init_fill_done && !full_reg)) begin
      // shift register of past M samples
      for (int i = M-1; i > 0; i--)
        prev_samples[i] <= prev_samples[i-1];
      prev_samples[0] <= data_in;

      // stash current window
      for (int m = 0; m <= M; m++)
        row_buf[m][wr_ptr] <= (m == 0) ? data_in : prev_samples[m-1];

      wr_ptr <= wr_ptr + 1;
    end else if (full) begin
      wr_ptr <= rd_ptr;
    end
  end

  //----------------------------------------------------------------------  
  // READ: one cycle per handshake, when not busy & not empty
  //----------------------------------------------------------------------  
  wire do_read = (!read_issued && !mac_busy && (count >= WIDTH));

  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      rd_ptr      <= 0;
      data_valid  <= 0;
      read_issued <= 0;
    end else begin
      if (do_read) begin
        // output window
        for (int i = 0; i <= M; i++)
          data_out[i] <= row_buf[i][rd_ptr];

        rd_ptr      <= rd_ptr + 1;
        data_valid  <= 1;
        read_issued <= 1;
      end else begin
        data_valid  <= 0;
        read_issued <= 0;
        if (rd_ptr == WIDTH)
          rd_ptr <= 0;
      end

      // clear issued flag when MAC goes busy again
      if (mac_busy)
        read_issued <= 0;
    end
  end

  //------------------------------------------------------------------------
  // COUNT: tracks writes minus reads
  //------------------------------------------------------------------------
  always_ff @(posedge clk or negedge rst_n) begin
    if (!rst_n) begin
      count <= 0;
    end else begin
      inc   = (valid_in && (count < WIDTH) && (!full && !init_fill_done_reg || init_fill_done && !full_reg)) ? 1 : 0;
      dec   = (do_read)                     ? 1 : 0;
      count <= count + inc - dec;
    end
  end

endmodule

  // ------------------------------------------------------------------------
  // Output assignment
  // ------------------------------------------------------------------------
  // always_comb begin
  //   if (data_valid && !mac_busy) begin
  //     data_out[0] = row0_buf[rd_ptr];
  //     data_out[1] = row1_buf[rd_ptr];
  //   end else begin
  //     data_out[0] = '0;
  //     data_out[1] = '0;
  //   end
  // end

  // // mux out current column of each row
  // always_comb begin
  //   if (data_valid) begin
  //     data_out[0] = row0[rd_ptr];
  //     data_out[1] = row1[rd_ptr];
  //   end else begin
  //     data_out[0] = '0;
  //     data_out[1] = '0;
  //   end
  // end
  //end
// always_ff @(posedge clk or negedge rst_n) begin : proc_
//     if(!rst_n) begin
//        data_out[0] <= '0;
//       data_out[1] <= '0;  
//     end else begin
//          if (data_valid) begin
//       data_out[0] <= row0[rd_ptr];
//       data_out[1] <= row1[rd_ptr];
//     end
// end
// end


// endmodule





//////////////////////////////////////works but no width/////////////////////////////
// module sample_buffer #(
//     parameter DATA_WIDTH = 16,  // bit-width of the input samples
//     parameter M = 1             // memory depth (number of past samples to store)
// )(
//     input  logic                    clk,
//     input  logic                    rst,     // asynchronous active-low reset
//     input  logic                    valid_in,  // new sample valid
//     input logic                     mac_busy,
//     input  logic signed [DATA_WIDTH-1:0] data_in,
//     output logic                    data_valid,  // high when data_out is valid for MAC
//     output logic signed [DATA_WIDTH-1:0] data_out [0:M]  // output buffer window: [0]=current, [1]=prev, etc.
// );

//     // Internal shift register buffer
//     logic signed [DATA_WIDTH-1:0] buffer [0:M];
//        logic                         valid_reg;

//     // Shift logic
//     always_ff @(posedge clk or negedge rst) begin
//         if (!rst) begin
//             for (int i = 0; i <= M; i++) begin
//                 buffer[i] <= '0;
//             end
//         end else if (valid_in && !mac_busy) begin
//             // Shift older samples
//             for (int i = M; i > 0; i--) begin
//                 buffer[i] <= buffer[i-1];
//             end
//             // Store new input
//             buffer[0] <= data_in;
//              valid_reg <= 1'b1;
//         end
//     end
//    // Data valid is asserted when buffer has at least one sample and MAC is not busy
//     assign data_valid = valid_reg && !mac_busy;
//     // Output assignment
//     always_comb begin
//         for (int i = 0; i <= M; i++) begin
//             if(!mac_busy)
//                 data_out[i] = buffer[i];

//         end
//     end

// endmodule
