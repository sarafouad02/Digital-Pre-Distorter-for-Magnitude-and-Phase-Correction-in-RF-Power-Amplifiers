// Dual-port coefficient ROM (read-only), now loading from .mem files
module coeff_rom #(
  parameter COEFF_WIDTH    = 16,
  parameter DEPTH          = 6,
  parameter FRAC_SZ        = 12,
  // filenames for readmemh
  parameter  REAL_MEM_FILE = "coeffs_real.mem",
  parameter  IMAG_MEM_FILE = "coeffs_imag.mem"
)(
  input  logic                          clk,
  input  logic [$clog2(DEPTH)-1:0]      addr_re,
  input  logic [$clog2(DEPTH)-1:0]      addr_im,

  input  logic                          enable,
  output logic signed [COEFF_WIDTH-1:0] data_re,
  output logic signed [COEFF_WIDTH-1:0] data_im
);

// 1) Add a Vivado attribute to force block‐RAM style
(* ram_style = "block" *)
logic signed [COEFF_WIDTH-1:0] rom_re [0:DEPTH-1];
(* ram_style = "block" *)
logic signed [COEFF_WIDTH-1:0] rom_im [0:DEPTH-1];
  // initialize at elaboration by loading hex words from files
  initial begin
    // If the files aren’t found, simulation will warn/error
    $display("Loading real coeffs from %s", REAL_MEM_FILE);
    $readmemh(REAL_MEM_FILE, rom_re);

    $display("Loading imag coeffs from %s", IMAG_MEM_FILE);
    $readmemh(IMAG_MEM_FILE, rom_im);
  end

  // synchronous read: register outputs when enabled
  always_ff @(posedge clk) begin
    if (enable) begin
      data_re <= rom_re[addr_re];
      data_im <= rom_im[addr_im];
    end
  end

endmodule
