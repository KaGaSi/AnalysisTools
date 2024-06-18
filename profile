--- Startup times for process: Primary/TUI ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.003  000.003: --- NVIM STARTING ---
000.276  000.273: event init
000.390  000.115: early init
000.458  000.068: locale set
000.521  000.063: init first window
000.890  000.368: inits 1
000.902  000.012: window checked
001.018  000.116: parsing arguments
001.417  000.025  000.025: require('vim.shared')
001.500  000.035  000.035: require('vim.inspect')
001.534  000.024  000.024: require('vim._options')
001.535  000.115  000.056: require('vim._editor')
001.535  000.172  000.031: require('vim._init_packages')
001.536  000.346: init lua interpreter
002.256  000.719: --- NVIM STARTED ---

--- Startup times for process: Embedded ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.136  000.134: event init
000.180  000.045: early init
000.210  000.030: locale set
000.238  000.028: init first window
000.387  000.150: inits 1
000.396  000.009: window checked
000.447  000.051: parsing arguments
000.771  000.024  000.024: require('vim.shared')
000.835  000.030  000.030: require('vim.inspect')
000.869  000.025  000.025: require('vim._options')
000.870  000.097  000.042: require('vim._editor')
000.871  000.137  000.017: require('vim._init_packages')
000.871  000.287: init lua interpreter
000.931  000.059: expanding arguments
000.965  000.034: inits 2
001.150  000.185: init highlight
001.151  000.001: waiting for UI
001.276  000.125: done waiting for UI
001.279  000.003: clear screen
001.324  000.003  000.003: require('vim.keymap')
001.816  000.535  000.532: require('vim._defaults')
001.817  000.003: init default mappings & autocommands
002.074  000.048  000.048: sourcing /usr/share/nvim/runtime/ftplugin.vim
002.118  000.015  000.015: sourcing /usr/share/nvim/runtime/indent.vim
002.212  000.062  000.062: sourcing /usr/share/nvim/archlinux.lua
002.215  000.080  000.018: sourcing /etc/xdg/nvim/sysinit.vim
003.196  000.757  000.757: require('KaGaSi.core.keymap')
004.006  000.808  000.808: require('KaGaSi.core.options')
004.072  000.060  000.060: require('KaGaSi.core.fold')
004.100  000.026  000.026: require('KaGaSi.core.filepos')
004.118  000.017  000.017: require('KaGaSi.core.pdf')
004.119  001.869  000.201: require('KaGaSi.core')
004.462  000.228  000.228: require('lazy')
004.477  000.007  000.007: require('ffi')
004.533  000.010  000.010: require('vim.fs')
004.598  000.063  000.063: require('vim.uri')
004.604  000.125  000.053: require('vim.loader')
004.959  000.343  000.343: require('lazy.stats')
005.061  000.083  000.083: require('lazy.core.util')
005.145  000.083  000.083: require('lazy.core.config')
005.301  000.055  000.055: require('lazy.core.handler')
005.397  000.095  000.095: require('lazy.core.plugin')
005.408  000.261  000.111: require('lazy.core.loader')
005.894  000.065  000.065: require('KaGaSi.plugins.alpha')
006.021  000.116  000.116: require('KaGaSi.plugins.auto-session')
006.104  000.071  000.071: require('KaGaSi.plugins.autopairs')
006.172  000.049  000.049: require('KaGaSi.plugins.bigfile')
006.270  000.084  000.084: require('KaGaSi.plugins.bufferline')
006.373  000.075  000.075: require('KaGaSi.plugins.colorscheme')
006.478  000.094  000.094: require('KaGaSi.plugins.comment')
006.524  000.036  000.036: require('KaGaSi.plugins.dressing')
006.560  000.030  000.030: require('KaGaSi.plugins.gitsigns')
006.589  000.024  000.024: require('KaGaSi.plugins.indent-blankline')
006.640  000.047  000.047: require('KaGaSi.plugins.leap')
006.678  000.033  000.033: require('KaGaSi.plugins.leap-spooky')
006.725  000.035  000.035: require('KaGaSi.plugins.lspconfig')
006.966  000.029  000.029: require('KaGaSi.plugins.ltex_extra')
007.012  000.043  000.043: require('KaGaSi.plugins.lualine')
007.094  000.073  000.073: require('KaGaSi.plugins.mason')
007.135  000.031  000.031: require('KaGaSi.plugins.numbertoggle')
007.171  000.031  000.031: require('KaGaSi.plugins.nvim-cmp')
007.234  000.046  000.046: require('KaGaSi.plugins.nvim-tree')
007.264  000.023  000.023: require('KaGaSi.plugins.surround')
007.313  000.027  000.027: require('KaGaSi.plugins.telescope')
007.360  000.024  000.024: require('KaGaSi.plugins.todo-comments')
007.386  000.020  000.020: require('KaGaSi.plugins.toggleterm')
007.456  000.052  000.052: require('KaGaSi.plugins.treesitter')
007.534  000.075  000.075: require('KaGaSi.plugins.trouble')
007.679  000.117  000.117: require('KaGaSi.plugins.vim-maximizer')
007.726  000.043  000.043: require('KaGaSi.plugins.vimtex')
008.786  001.023  001.023: require('vim.filetype')
008.789  001.056  000.033: require('KaGaSi.plugins.vimwiki')
008.840  000.046  000.046: require('KaGaSi.plugins.which_key')
009.151  000.031  000.031: require('lazy.core.handler.event')
009.153  000.072  000.041: require('lazy.core.handler.ft')
009.179  000.023  000.023: require('lazy.core.handler.cmd')
009.205  000.025  000.025: require('lazy.core.handler.keys')
009.441  000.021  000.021: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/cls.vim
009.471  000.011  000.011: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tex.vim
009.494  000.008  000.008: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tikz.vim
009.899  000.120  000.120: sourcing /usr/share/nvim/runtime/filetype.lua
010.123  000.120  000.120: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/scroll.vim
010.170  000.026  000.026: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/themer.lua
011.239  000.322  000.322: require('themer')
011.281  000.039  000.039: require('themer.config')
011.318  000.030  000.030: require('themer.modules.core')
011.375  000.055  000.055: require('themer.modules.core.api')
011.443  000.066  000.066: require('themer.modules.themes.sonokai_deep')
011.477  000.032  000.032: require('themer.modules.core.utils')
011.576  000.098  000.098: require('themer.modules.core.mapper')
013.017  002.140  001.497: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/colors/themer_sonokai_deep.lua
014.209  000.082  000.082: require('mason-core.functional')
014.252  000.041  000.041: require('mason-core.path')
014.289  000.036  000.036: require('mason.settings')
014.300  000.259  000.101: require('mason-core.log')
014.302  000.322  000.062: require('mason-core.EventEmitter')
014.342  000.039  000.039: require('mason-core.optional')
014.454  000.068  000.068: require('mason-core.async')
014.501  000.044  000.044: require('mason-core.async.uv')
014.507  000.164  000.051: require('mason-core.fs')
014.554  000.046  000.046: require('mason-registry.sources')
014.633  000.036  000.036: require('mason-core.functional.data')
014.677  000.043  000.043: require('mason-core.functional.function')
014.691  000.135  000.056: require('mason-core.functional.list')
014.701  000.955  000.249: require('mason-registry')
014.705  001.318  000.363: require('mason-tool-installer')
014.715  001.372  000.054: sourcing /home/gary/.local/share/nvim/lazy/mason-tool-installer.nvim/plugin/mason-tool-installer.lua
014.891  000.037  000.037: require('mason-core.functional.relation')
014.931  000.036  000.036: require('mason-core.functional.logic')
014.974  000.188  000.115: require('mason-core.platform')
014.977  000.230  000.041: require('mason')
015.057  000.033  000.033: require('mason-lspconfig.settings')
015.060  000.082  000.049: require('mason-lspconfig')
015.285  000.048  000.048: require('mason-core.functional.string')
015.304  000.148  000.100: require('mason.api.command')
015.422  000.035  000.035: require('mason-lspconfig.notify')
015.424  000.081  000.046: require('mason-lspconfig.lspconfig_hook')
016.948  000.031  000.031: sourcing /home/gary/.local/share/nvim/lazy/plenary.nvim/plugin/plenary.vim
017.259  000.035  000.035: require('plenary.bit')
017.290  000.028  000.028: require('plenary.functional')
017.314  000.291  000.228: require('plenary.path')
017.332  000.008  000.008: require('vim.F')
017.354  000.853  000.522: require('plenary.log')
017.387  001.029  000.176: require('lsp-file-operations.log')
017.388  001.341  000.312: require('lsp-file-operations')
017.658  000.021  000.021: sourcing /home/gary/.local/share/nvim/lazy/nvim-web-devicons/plugin/nvim-web-devicons.vim
018.086  000.023  000.023: require('nvim-tree.notify')
018.089  000.050  000.027: require('nvim-tree.events')
018.176  000.024  000.024: require('nvim-tree.log')
018.244  000.032  000.032: require('nvim-tree.iterators.node-iterator')
018.263  000.086  000.053: require('nvim-tree.utils')
018.291  000.027  000.027: require('nvim-tree.git.utils')
018.319  000.027  000.027: require('nvim-tree.git.runner')
018.347  000.027  000.027: require('nvim-tree.watcher')
018.373  000.024  000.024: require('nvim-tree.explorer.node')
018.375  000.260  000.045: require('nvim-tree.git')
018.398  000.022  000.022: require('nvim-tree.explorer.watch')
018.471  000.022  000.022: require('nvim-tree.explorer.node-builders')
018.495  000.023  000.023: require('nvim-tree.explorer.sorters')
018.543  000.022  000.022: require('nvim-tree.marks')
018.545  000.049  000.027: require('nvim-tree.explorer.filters')
018.612  000.041  000.041: require('nvim-tree.view')
018.615  000.069  000.028: require('nvim-tree.live-filter')
018.616  000.192  000.028: require('nvim-tree.explorer.explore')
018.652  000.035  000.035: require('nvim-tree.explorer.reload')
018.654  000.564  000.055: require('nvim-tree.explorer')
018.656  000.638  000.023: require('nvim-tree.core')
018.680  000.023  000.023: require('nvim-tree.renderer.components.padding')
018.704  000.023  000.023: require('nvim-tree.renderer.components.icons')
018.729  000.024  000.024: require('nvim-tree.renderer.components.full-name')
018.798  000.018  000.018: require('nvim-tree.enum')
018.821  000.022  000.022: require('nvim-tree.renderer.decorator')
018.822  000.061  000.022: require('nvim-tree.renderer.decorator.bookmarks')
018.842  000.019  000.019: require('nvim-tree.renderer.decorator.copied')
018.862  000.019  000.019: require('nvim-tree.renderer.decorator.cut')
018.926  000.033  000.033: require('nvim-tree.diagnostics')
019.083  000.155  000.155: require('vim.diagnostic')
019.086  000.223  000.035: require('nvim-tree.renderer.decorator.diagnostics')
019.142  000.056  000.056: require('nvim-tree.renderer.decorator.git')
019.338  000.164  000.164: require('nvim-tree.buffers')
019.340  000.197  000.033: require('nvim-tree.renderer.decorator.modified')
019.367  000.026  000.026: require('nvim-tree.renderer.decorator.opened')
019.379  000.649  000.049: require('nvim-tree.renderer.builder')
019.382  001.404  000.047: require('nvim-tree.renderer')
019.386  001.441  000.037: require('nvim-tree.lib')
019.450  000.063  000.063: require('nvim-tree.appearance')
019.803  000.036  000.036: require('nvim-tree.actions.finders.find-file')
019.846  000.041  000.041: require('nvim-tree.actions.finders.search-node')
019.847  000.112  000.035: require('nvim-tree.actions.finders')
019.984  000.052  000.052: require('nvim-tree.actions.reloaders')
019.990  000.111  000.059: require('nvim-tree.actions.fs.copy-paste')
020.033  000.042  000.042: require('nvim-tree.actions.fs.create-file')
020.077  000.043  000.043: require('nvim-tree.actions.fs.remove-file')
020.122  000.044  000.044: require('nvim-tree.actions.fs.rename-file')
020.168  000.045  000.045: require('nvim-tree.actions.fs.trash')
020.170  000.323  000.038: require('nvim-tree.actions.fs')
020.244  000.039  000.039: require('nvim-tree.actions.moves.item')
020.287  000.042  000.042: require('nvim-tree.actions.moves.parent')
020.319  000.031  000.031: require('nvim-tree.actions.moves.sibling')
020.320  000.149  000.036: require('nvim-tree.actions.moves')
020.393  000.038  000.038: require('nvim-tree.actions.node.file-popup')
020.466  000.072  000.072: require('nvim-tree.actions.node.open-file')
020.501  000.033  000.033: require('nvim-tree.actions.node.run-command')
020.535  000.034  000.034: require('nvim-tree.actions.node.system-open')
020.537  000.216  000.039: require('nvim-tree.actions.node')
020.622  000.054  000.054: require('nvim-tree.actions.root.change-dir')
020.656  000.033  000.033: require('nvim-tree.actions.root.dir-up')
020.657  000.119  000.033: require('nvim-tree.actions.root')
020.742  000.038  000.038: require('nvim-tree.actions.tree.find-file')
020.810  000.033  000.033: require('nvim-tree.actions.tree.modifiers.collapse-all')
020.852  000.041  000.041: require('nvim-tree.actions.tree.modifiers.expand-all')
020.889  000.035  000.035: require('nvim-tree.actions.tree.modifiers.toggles')
020.890  000.147  000.037: require('nvim-tree.actions.tree.modifiers')
020.922  000.031  000.031: require('nvim-tree.actions.tree.open')
020.955  000.033  000.033: require('nvim-tree.actions.tree.toggle')
020.957  000.298  000.050: require('nvim-tree.actions.tree')
020.958  001.275  000.058: require('nvim-tree.actions')
021.019  000.061  000.061: require('nvim-tree.appearance.diagnostics')
021.205  000.116  000.116: require('nvim-tree.keymap')
021.231  000.211  000.095: require('nvim-tree.help')
021.300  000.068  000.068: require('nvim-tree.marks.navigation')
021.387  000.087  000.087: require('nvim-tree.marks.bulk-delete')
021.486  000.097  000.097: require('nvim-tree.marks.bulk-trash')
021.599  000.112  000.112: require('nvim-tree.marks.bulk-move')
021.649  002.160  000.251: require('nvim-tree.api')
021.659  002.209  000.048: require('nvim-tree.commands')
021.720  000.059  000.059: require('nvim-tree.legacy')
021.735  004.037  000.265: require('nvim-tree')
022.633  000.260  000.260: require('nvim-web-devicons.icons-default')
022.758  000.447  000.187: require('nvim-web-devicons')
024.483  007.091  002.586: require('nvim-tree.api')
025.496  000.131  000.131: require('neodev')
025.523  000.025  000.025: require('neodev.config')
025.576  000.023  000.023: require('neodev.util')
025.578  000.050  000.027: require('neodev.lsp')
025.762  000.059  000.059: require('vim.lsp.log')
026.047  000.283  000.283: require('vim.lsp.protocol')
026.344  000.137  000.137: require('vim.lsp._snippet_grammar')
026.383  000.037  000.037: require('vim.highlight')
026.398  000.348  000.174: require('vim.lsp.util')
026.498  000.056  000.056: require('vim.lsp.sync')
026.506  000.105  000.049: require('vim.lsp._changetracking')
026.577  000.070  000.070: require('vim.lsp.rpc')
026.620  001.006  000.140: require('vim.lsp')
026.679  001.100  000.094: require('lspconfig.util')
026.995  000.067  000.067: sourcing /home/gary/.local/share/nvim/lazy/nvim-lspconfig/plugin/lspconfig.lua
027.137  000.035  000.035: require('lspconfig.async')
027.139  000.097  000.062: require('lspconfig.configs')
027.141  000.133  000.036: require('lspconfig')
027.238  000.055  000.055: require('cmp_nvim_lsp.source')
027.240  000.098  000.043: require('cmp_nvim_lsp')
027.385  000.053  000.053: require('mason-core.functional.table')
027.415  000.151  000.098: require('mason-lspconfig.mappings.server')
027.532  000.053  000.053: require('lspconfig.server_configurations.pyright')
027.853  000.073  000.073: require('lspconfig.manager')
027.904  000.046  000.046: require('lspconfig.server_configurations.cssls')
028.110  000.052  000.052: require('lspconfig.server_configurations.html')
028.308  000.101  000.101: require('lspconfig.server_configurations.ltex')
028.460  000.043  000.043: require('lspconfig.server_configurations.bashls')
028.576  000.046  000.046: require('lspconfig.server_configurations.vimls')
028.685  000.048  000.048: require('lspconfig.server_configurations.lua_ls')
028.793  000.046  000.046: require('lspconfig.server_configurations.clangd')
028.874  013.449  002.753: require('lspconfig.util')
028.919  000.042  000.042: require('mason-lspconfig.server_config_extensions')
028.975  000.055  000.055: require('lspconfig.server_configurations.omnisharp')
029.027  000.040  000.040: require('mason-lspconfig.ensure_installed')
029.149  000.067  000.067: require('mason-core.result')
029.330  000.087  000.087: require('mason-core.process')
029.464  000.133  000.133: require('mason-core.spawn')
029.475  000.281  000.061: require('mason-core.managers.powershell')
029.510  000.033  000.033: require('mason.version')
029.513  000.362  000.048: require('mason-core.fetch')
029.562  000.046  000.046: require('mason-core.providers')
029.828  000.146  000.146: require('mason-core.purl')
029.853  000.239  000.094: require('mason-core.package')
030.115  000.117  000.117: require('mason-core.installer.registry.expr')
030.132  000.202  000.085: require('mason-core.installer.registry.link')
031.083  000.260  000.260: require('mason-core.receipt')
031.109  000.630  000.370: require('mason-core.installer.context')
031.202  000.091  000.091: require('mason-core.async.control')
031.281  000.077  000.077: require('mason-core.installer.linker')
031.288  000.971  000.173: require('mason-core.installer')
031.303  001.091  000.120: require('mason-core.installer.managers.std')
031.305  001.171  000.081: require('mason-core.installer.registry.schemas')
031.391  000.039  000.039: require('mason-core.installer.registry.util')
031.400  001.545  000.133: require('mason-core.installer.registry')
031.401  001.838  000.054: require('mason-registry.sources.util')
031.406  002.373  000.060: require('mason-registry.sources.github')
034.811  000.025  000.025: require('mason-core.functional.number')
034.824  000.078  000.053: require('mason-lspconfig.api.command')
035.520  000.038  000.038: require('lualine_require')
035.685  000.373  000.335: require('lualine')
035.713  000.026  000.026: require('lazy.status')
036.700  000.029  000.029: require('lualine.utils.mode')
039.642  000.020  000.020: sourcing /home/gary/.local/share/nvim/lazy/vim-numbertoggle/plugin/number_toggle.vim
040.068  000.236  000.236: sourcing /home/gary/.local/share/nvim/lazy/leap.nvim/plugin/init.lua
040.438  000.186  000.186: require('leap-spooky')
041.417  000.040  000.040: require('auto-session.logger')
041.423  000.081  000.041: require('auto-session.lib')
041.453  000.029  000.029: require('auto-session.autocmds')
041.495  000.367  000.256: require('auto-session')
042.060  000.014  000.014: sourcing /home/gary/.local/share/nvim/lazy/todo-comments.nvim/plugin/todo.vim
042.306  000.025  000.025: require('todo-comments.util')
042.313  000.075  000.049: require('todo-comments.config')
042.376  000.039  000.039: require('todo-comments.highlight')
042.377  000.064  000.024: require('todo-comments.jump')
042.379  000.307  000.169: require('todo-comments')
042.587  000.124  000.124: sourcing /home/gary/.local/share/nvim/lazy/telescope.nvim/plugin/telescope.lua
042.679  000.038  000.038: require('telescope._extensions')
042.681  000.081  000.043: require('telescope')
042.954  000.052  000.052: require('plenary.strings')
042.990  000.034  000.034: require('telescope.deprecated')
043.134  000.083  000.083: require('telescope.log')
043.290  000.026  000.026: require('plenary.compat')
043.296  000.089  000.063: require('plenary.job')
043.331  000.034  000.034: require('telescope.state')
043.350  000.215  000.093: require('telescope.utils')
043.355  000.364  000.066: require('telescope.sorters')
044.292  001.470  001.019: require('telescope.config')
044.473  000.111  000.111: require('plenary.window.border')
044.505  000.031  000.031: require('plenary.window')
044.532  000.026  000.026: require('plenary.popup.utils')
044.534  000.240  000.072: require('plenary.popup')
044.580  000.044  000.044: require('telescope.pickers.scroller')
044.619  000.038  000.038: require('telescope.actions.state')
044.663  000.043  000.043: require('telescope.actions.utils')
044.762  000.052  000.052: require('telescope.actions.mt')
044.775  000.111  000.060: require('telescope.actions.set')
044.862  000.046  000.046: require('telescope.config.resolve')
044.864  000.088  000.042: require('telescope.pickers.entry_display')
044.903  000.038  000.038: require('telescope.from_entry')
045.384  002.702  000.630: require('telescope.actions')
047.035  000.299  000.299: require('fzf_lib')
047.047  000.489  000.190: require('telescope._extensions.fzf')
047.135  005.627  001.911: require('telescope')
047.523  000.228  000.228: require('auto-session.session-lens.library')
047.711  000.185  000.185: require('auto-session.session-lens.actions')
047.725  000.588  000.175: require('auto-session.session-lens')
048.336  000.080  000.080: require('bufferline.lazy')
048.493  000.147  000.147: require('bufferline.commands')
048.748  000.253  000.253: require('bufferline.config')
048.763  000.840  000.361: require('bufferline')
049.084  000.138  000.138: require('bufferline.utils')
049.086  000.314  000.176: require('bufferline.groups')
049.176  000.071  000.071: require('bufferline.constants')
049.296  000.116  000.116: require('bufferline.colors')
049.769  000.207  000.207: require('bufferline.highlights')
051.647  000.135  000.135: require('vim.version')
053.004  001.566  001.431: require('bufferline.hover')
053.510  000.339  000.339: require('bufferline.ui')
054.121  000.104  000.104: require('toggleterm.lazy')
054.166  000.042  000.042: require('toggleterm.constants')
054.361  000.193  000.193: require('toggleterm.terminal')
054.371  000.648  000.309: require('toggleterm')
054.472  000.043  000.043: require('toggleterm.colors')
054.717  000.244  000.244: require('toggleterm.utils')
054.729  000.358  000.071: require('toggleterm.config')
056.237  000.253  000.253: require('toggleterm.commandline')
056.739  000.133  000.133: sourcing /usr/share/nvim/runtime/plugin/editorconfig.lua
056.862  000.090  000.090: sourcing /usr/share/nvim/runtime/plugin/gzip.vim
057.198  000.293  000.293: sourcing /usr/share/nvim/runtime/plugin/man.lua
057.591  000.115  000.115: sourcing /usr/share/nvim/runtime/pack/dist/opt/matchit/plugin/matchit.vim
057.602  000.380  000.266: sourcing /usr/share/nvim/runtime/plugin/matchit.vim
057.697  000.064  000.064: sourcing /usr/share/nvim/runtime/plugin/matchparen.vim
057.737  000.008  000.008: sourcing /usr/share/nvim/runtime/plugin/netrwPlugin.vim
057.912  000.162  000.162: sourcing /usr/share/nvim/runtime/plugin/osc52.lua
058.040  000.097  000.097: sourcing /usr/share/nvim/runtime/plugin/rplugin.vim
058.096  000.029  000.029: sourcing /usr/share/nvim/runtime/plugin/shada.vim
058.143  000.010  000.010: sourcing /usr/share/nvim/runtime/plugin/spellfile.vim
058.226  000.045  000.045: sourcing /usr/share/nvim/runtime/plugin/tarPlugin.vim
058.355  000.096  000.096: sourcing /usr/share/nvim/runtime/plugin/tohtml.lua
058.388  000.012  000.012: sourcing /usr/share/nvim/runtime/plugin/tutor.vim
058.485  000.062  000.062: sourcing /usr/share/nvim/runtime/plugin/zipPlugin.vim
059.231  000.085  000.085: require('bigfile.features')
059.234  000.449  000.364: require('bigfile')
059.247  000.530  000.081: sourcing /home/gary/.local/share/nvim/lazy/bigfile.nvim/after/plugin/bigfile.lua
059.312  000.025  000.025: sourcing /home/gary/.local/share/nvim/lazy/cmp-nvim-lsp/after/plugin/cmp_nvim_lsp.lua
059.340  055.221  016.884: require('KaGaSi.lazy')
059.341  057.109  000.018: sourcing /home/gary/.config/nvim/init.lua
059.346  000.276: sourcing vimrc file(s)
059.485  000.021  000.021: sourcing /usr/share/nvim/runtime/filetype.lua
059.634  000.035  000.035: sourcing /usr/share/nvim/runtime/syntax/synload.vim
060.226  000.704  000.669: sourcing /usr/share/nvim/runtime/syntax/syntax.vim
060.235  000.165: inits 3
061.565  001.329: reading ShaDa
061.583  000.019: making windows
062.250  000.041  000.041: require('ibl.utils')
062.256  000.090  000.050: require('ibl.config')
062.329  000.043  000.043: require('ibl.indent')
062.334  000.078  000.034: require('ibl.hooks')
062.337  000.199  000.031: require('ibl.highlights')
062.374  000.035  000.035: require('ibl.autocmds')
062.497  000.121  000.121: require('ibl.inlay_hints')
062.527  000.028  000.028: require('ibl.virt_text')
062.646  000.096  000.096: require('ibl.scope_languages')
062.648  000.120  000.024: require('ibl.scope')
062.651  000.787  000.284: require('ibl')
062.663  000.827  000.040: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/plugin/commands.lua
062.834  000.062  000.062: require('vim.iter')
062.922  000.076  000.076: require('vim.lsp.handlers')
063.217  000.035  000.035: require('vim.treesitter.language')
063.239  000.020  000.020: require('vim.func')
063.265  000.024  000.024: require('vim.func._memoize')
063.274  000.146  000.067: require('vim.treesitter.query')
063.301  000.026  000.026: require('vim.treesitter._range')
063.309  000.226  000.055: require('vim.treesitter.languagetree')
063.313  000.275  000.049: require('vim.treesitter')
064.529  000.038  000.038: require('nvim-treesitter.utils')
065.031  000.500  000.500: require('nvim-treesitter.parsers')
065.161  000.022  000.022: require('nvim-treesitter.compat')
065.218  000.031  000.031: require('nvim-treesitter.ts_utils')
065.221  000.058  000.027: require('nvim-treesitter.tsrange')
065.242  000.020  000.020: require('nvim-treesitter.caching')
065.247  000.136  000.035: require('nvim-treesitter.query')
065.253  000.179  000.044: require('nvim-treesitter.configs')
065.255  000.222  000.042: require('nvim-treesitter.info')
065.308  000.037  000.037: require('nvim-treesitter.shell_command_selectors')
065.323  000.892  000.095: require('nvim-treesitter.install')
065.347  000.023  000.023: require('nvim-treesitter.statusline')
065.377  000.029  000.029: require('nvim-treesitter.query_predicates')
065.378  001.213  000.269: require('nvim-treesitter')
065.613  001.489  000.275: sourcing /home/gary/.local/share/nvim/lazy/nvim-treesitter/plugin/nvim-treesitter.lua
067.052  000.033  000.033: require('nvim-treesitter.locals')
067.056  000.073  000.041: require('nvim-treesitter.incremental_selection')
067.288  000.030  000.030: require('nvim-treesitter.highlight')
067.417  000.038  000.038: require('nvim-treesitter.indent')
067.726  000.070  000.070: sourcing /home/gary/.local/share/nvim/lazy/nvim-ts-context-commentstring/plugin/ts_context_commentstring.lua
068.018  000.024  000.024: require('Comment.config')
068.096  000.043  000.043: require('Comment.ft')
068.098  000.079  000.036: require('Comment.utils')
068.124  000.025  000.025: require('Comment.opfunc')
068.149  000.024  000.024: require('Comment.extra')
068.151  000.357  000.205: require('Comment.api')
068.243  000.472  000.115: sourcing /home/gary/.local/share/nvim/lazy/Comment.nvim/plugin/Comment.lua
068.286  000.027  000.027: require('Comment')
068.312  000.024  000.024: require('ts_context_commentstring.integrations.comment_nvim')
068.781  000.027  000.027: require('gitsigns.async')
068.807  000.023  000.023: require('gitsigns.debug.log')
068.860  000.053  000.053: require('gitsigns.config')
068.863  000.338  000.235: require('gitsigns')
068.938  000.042  000.042: require('gitsigns.highlight')
069.692  000.105  000.105: require('gitsigns.util')
069.749  000.052  000.052: require('gitsigns.system')
069.800  000.048  000.048: require('gitsigns.message')
069.857  000.055  000.055: require('gitsigns.git.version')
069.872  000.426  000.167: require('gitsigns.git')
070.037  000.064  000.064: require('gitsigns.cache')
070.096  000.056  000.056: require('gitsigns.signs')
070.147  000.049  000.049: require('gitsigns.status')
070.195  000.046  000.046: require('gitsigns.debounce')
070.242  000.044  000.044: require('gitsigns.diff')
070.323  000.080  000.080: require('gitsigns.hunks')
070.337  000.463  000.123: require('gitsigns.manager')
070.345  001.011  000.123: require('gitsigns.attach')
070.400  000.042  000.042: require('gitsigns.current_line_blame')
070.477  000.043  000.043: require('vim._system')
072.827  000.073  000.073: require('nvim-surround.input')
072.884  000.052  000.052: require('nvim-surround.functional')
072.901  000.529  000.404: require('nvim-surround.config')
072.916  000.593  000.065: require('nvim-surround.buffer')
072.971  000.053  000.053: require('nvim-surround.cache')
073.047  000.074  000.074: require('nvim-surround.utils')
073.053  001.071  000.351: require('nvim-surround')
076.146  000.145  000.145: sourcing /usr/share/nvim/runtime/autoload/provider/clipboard.vim
077.364  000.264  000.264: require('vim.filetype.detect')
077.975  000.053  000.053: sourcing /usr/share/nvim/runtime/ftplugin/c.vim
078.037  000.048  000.048: sourcing /usr/share/nvim/runtime/ftplugin/c.lua
078.311  000.044  000.044: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/ftplugin/c.lua
078.596  000.012  000.012: sourcing /usr/share/nvim/runtime/indent/c.vim
079.065  000.121  000.121: require('vim.lsp.client')
079.427  000.089  000.089: require('vim.glob')
079.434  000.181  000.092: require('vim.lsp._dynamic')
085.226  000.579  000.579: sourcing /usr/share/nvim/runtime/syntax/c.vim
126.835  000.162  000.162: require('vim.treesitter.highlighter')
127.705  000.240  000.240: require('ts_context_commentstring.utils')
127.903  000.196  000.196: require('ts_context_commentstring.config')
127.913  000.528  000.092: require('ts_context_commentstring.internal')
128.298  000.314  000.314: require('editorconfig')
128.506  058.461: opening buffers
133.621  000.963  000.963: require('bufferline.state')
133.642  004.172: BufEnter autocommands
138.679  000.027  000.027: sourcing /usr/share/nvim/runtime/ftplugin/conf.vim
139.950  000.039  000.039: sourcing /usr/share/nvim/runtime/syntax/conf.vim
152.754  019.046: editing files in windows
157.526  004.773: executing command arguments
158.166  000.436  000.436: require('alpha')
158.231  000.063  000.063: require('alpha.themes.dashboard')
169.450  011.425: VimEnter autocommands
169.518  000.068: UIEnter autocommands
169.531  000.014: diff scrollbinding
169.533  000.002: before starting main loop
170.176  000.069  000.069: require('bufferline.tabpages')
170.234  000.052  000.052: require('bufferline.models')
170.266  000.031  000.031: require('bufferline.pick')
170.295  000.027  000.027: require('bufferline.duplicates')
170.406  000.064  000.064: require('bufferline.diagnostics')
170.642  000.067  000.067: require('bufferline.numbers')
170.763  000.036  000.036: require('bufferline.sorters')
170.834  000.048  000.048: require('bufferline.offset')
170.889  000.052  000.052: require('bufferline.custom_area')
209.111  039.132: first screen update
209.119  000.009: --- NVIM STARTED ---

--- Startup times for process: Primary/TUI ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.003  000.003: --- NVIM STARTING ---
000.263  000.260: event init
000.391  000.128: early init
000.500  000.109: locale set
000.595  000.095: init first window
001.069  000.474: inits 1
001.084  000.015: window checked
001.190  000.106: parsing arguments
001.648  000.032  000.032: require('vim.shared')
001.719  000.038  000.038: require('vim.inspect')
001.763  000.033  000.033: require('vim._options')
001.764  000.114  000.043: require('vim._editor')
001.765  000.177  000.031: require('vim._init_packages')
001.769  000.403: init lua interpreter
002.275  000.506: --- NVIM STARTED ---

--- Startup times for process: Embedded ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.115  000.114: event init
000.161  000.045: early init
000.186  000.025: locale set
000.214  000.028: init first window
000.362  000.147: inits 1
000.371  000.010: window checked
000.425  000.053: parsing arguments
000.737  000.024  000.024: require('vim.shared')
000.794  000.027  000.027: require('vim.inspect')
000.829  000.027  000.027: require('vim._options')
000.830  000.091  000.038: require('vim._editor')
000.831  000.131  000.015: require('vim._init_packages')
000.832  000.276: init lua interpreter
000.894  000.063: expanding arguments
000.931  000.036: inits 2
001.085  000.154: init highlight
001.086  000.001: waiting for UI
001.230  000.144: done waiting for UI
001.232  000.002: clear screen
001.284  000.003  000.003: require('vim.keymap')
001.657  000.423  000.420: require('vim._defaults')
001.659  000.003: init default mappings & autocommands
001.919  000.048  000.048: sourcing /usr/share/nvim/runtime/ftplugin.vim
001.964  000.020  000.020: sourcing /usr/share/nvim/runtime/indent.vim
002.046  000.047  000.047: sourcing /usr/share/nvim/archlinux.lua
002.049  000.067  000.020: sourcing /etc/xdg/nvim/sysinit.vim
002.655  000.376  000.376: require('KaGaSi.core.keymap')
003.465  000.808  000.808: require('KaGaSi.core.options')
003.534  000.067  000.067: require('KaGaSi.core.fold')
003.572  000.036  000.036: require('KaGaSi.core.filepos')
003.599  000.026  000.026: require('KaGaSi.core.pdf')
003.600  001.525  000.212: require('KaGaSi.core')
003.914  000.239  000.239: require('lazy')
003.937  000.011  000.011: require('ffi')
003.984  000.020  000.020: require('vim.fs')
004.100  000.111  000.111: require('vim.uri')
004.106  000.167  000.037: require('vim.loader')
004.276  000.157  000.157: require('lazy.stats')
004.347  000.050  000.050: require('lazy.core.util')
004.415  000.067  000.067: require('lazy.core.config')
004.562  000.070  000.070: require('lazy.core.handler')
004.696  000.133  000.133: require('lazy.core.plugin')
004.708  000.292  000.088: require('lazy.core.loader')
005.154  000.128  000.128: require('KaGaSi.plugins.alpha')
005.301  000.104  000.104: require('KaGaSi.plugins.auto-session')
005.378  000.066  000.066: require('KaGaSi.plugins.autopairs')
005.447  000.059  000.059: require('KaGaSi.plugins.bigfile')
005.495  000.041  000.041: require('KaGaSi.plugins.bufferline')
005.548  000.045  000.045: require('KaGaSi.plugins.colorscheme')
005.582  000.029  000.029: require('KaGaSi.plugins.comment')
005.635  000.036  000.036: require('KaGaSi.plugins.dressing')
005.766  000.091  000.091: require('KaGaSi.plugins.gitsigns')
005.900  000.089  000.089: require('KaGaSi.plugins.indent-blankline')
006.171  000.236  000.236: require('KaGaSi.plugins.leap')
006.292  000.107  000.107: require('KaGaSi.plugins.leap-spooky')
006.412  000.089  000.089: require('KaGaSi.plugins.lspconfig')
006.704  000.059  000.059: require('KaGaSi.plugins.ltex_extra')
006.762  000.052  000.052: require('KaGaSi.plugins.lualine')
006.823  000.047  000.047: require('KaGaSi.plugins.mason')
006.879  000.045  000.045: require('KaGaSi.plugins.numbertoggle')
006.953  000.065  000.065: require('KaGaSi.plugins.nvim-cmp')
007.015  000.040  000.040: require('KaGaSi.plugins.nvim-tree')
007.065  000.030  000.030: require('KaGaSi.plugins.surround')
007.108  000.031  000.031: require('KaGaSi.plugins.telescope')
007.159  000.040  000.040: require('KaGaSi.plugins.todo-comments')
007.281  000.031  000.031: require('KaGaSi.plugins.toggleterm')
007.350  000.062  000.062: require('KaGaSi.plugins.treesitter')
007.430  000.036  000.036: require('KaGaSi.plugins.trouble')
007.464  000.026  000.026: require('KaGaSi.plugins.vim-maximizer')
007.498  000.029  000.029: require('KaGaSi.plugins.vimtex')
008.484  000.952  000.952: require('vim.filetype')
008.487  000.985  000.034: require('KaGaSi.plugins.vimwiki')
008.532  000.031  000.031: require('KaGaSi.plugins.which_key')
008.777  000.037  000.037: require('lazy.core.handler.event')
008.779  000.070  000.033: require('lazy.core.handler.ft')
008.813  000.032  000.032: require('lazy.core.handler.keys')
008.841  000.026  000.026: require('lazy.core.handler.cmd')
009.325  000.016  000.016: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/cls.vim
009.352  000.010  000.010: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tex.vim
009.376  000.008  000.008: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tikz.vim
009.591  000.136  000.136: sourcing /usr/share/nvim/runtime/filetype.lua
009.802  000.105  000.105: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/scroll.vim
009.854  000.031  000.031: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/themer.lua
010.952  000.353  000.353: require('themer')
011.026  000.072  000.072: require('themer.config')
011.119  000.040  000.040: require('themer.modules.core')
011.183  000.062  000.062: require('themer.modules.core.api')
011.252  000.067  000.067: require('themer.modules.themes.sonokai_deep')
011.291  000.037  000.037: require('themer.modules.core.utils')
011.406  000.114  000.114: require('themer.modules.core.mapper')
012.867  002.297  001.552: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/colors/themer_sonokai_deep.lua
013.077  000.028  000.028: sourcing /home/gary/.local/share/nvim/lazy/nvim-web-devicons/plugin/nvim-web-devicons.vim
014.249  000.070  000.070: require('nvim-tree.notify')
014.259  000.150  000.081: require('nvim-tree.events')
014.469  000.068  000.068: require('nvim-tree.log')
014.777  000.207  000.207: require('nvim-tree.iterators.node-iterator')
014.820  000.349  000.142: require('nvim-tree.utils')
014.936  000.115  000.115: require('nvim-tree.git.utils')
015.004  000.067  000.067: require('nvim-tree.git.runner')
015.036  000.031  000.031: require('nvim-tree.watcher')
015.062  000.025  000.025: require('nvim-tree.explorer.node')
015.065  000.744  000.090: require('nvim-tree.git')
015.090  000.024  000.024: require('nvim-tree.explorer.watch')
015.136  000.023  000.023: require('nvim-tree.explorer.node-builders')
015.162  000.025  000.025: require('nvim-tree.explorer.sorters')
015.234  000.026  000.026: require('nvim-tree.marks')
015.236  000.073  000.047: require('nvim-tree.explorer.filters')
015.302  000.039  000.039: require('nvim-tree.view')
015.304  000.067  000.029: require('nvim-tree.live-filter')
015.306  000.215  000.026: require('nvim-tree.explorer.explore')
015.331  000.025  000.025: require('nvim-tree.explorer.reload')
015.333  001.073  000.064: require('nvim-tree.explorer')
015.335  001.277  000.054: require('nvim-tree.core')
015.358  000.022  000.022: require('nvim-tree.renderer.components.padding')
015.382  000.023  000.023: require('nvim-tree.renderer.components.icons')
015.406  000.023  000.023: require('nvim-tree.renderer.components.full-name')
015.488  000.019  000.019: require('nvim-tree.enum')
015.522  000.032  000.032: require('nvim-tree.renderer.decorator')
015.528  000.085  000.034: require('nvim-tree.renderer.decorator.bookmarks')
015.552  000.023  000.023: require('nvim-tree.renderer.decorator.copied')
015.582  000.029  000.029: require('nvim-tree.renderer.decorator.cut')
015.649  000.035  000.035: require('nvim-tree.diagnostics')
015.739  000.003  000.003: require('vim.F')
015.760  000.109  000.106: require('vim.diagnostic')
015.762  000.180  000.035: require('nvim-tree.renderer.decorator.diagnostics')
015.817  000.054  000.054: require('nvim-tree.renderer.decorator.git')
015.866  000.020  000.020: require('nvim-tree.buffers')
015.868  000.050  000.030: require('nvim-tree.renderer.decorator.modified')
015.899  000.031  000.031: require('nvim-tree.renderer.decorator.opened')
015.905  000.498  000.048: require('nvim-tree.renderer.builder')
015.907  001.916  000.072: require('nvim-tree.renderer')
015.910  002.118  000.202: require('nvim-tree.lib')
015.980  000.068  000.068: require('nvim-tree.appearance')
016.160  000.038  000.038: require('nvim-tree.actions.finders.find-file')
016.218  000.056  000.056: require('nvim-tree.actions.finders.search-node')
016.242  000.143  000.050: require('nvim-tree.actions.finders')
016.391  000.029  000.029: require('nvim-tree.actions.reloaders')
016.397  000.079  000.050: require('nvim-tree.actions.fs.copy-paste')
016.449  000.051  000.051: require('nvim-tree.actions.fs.create-file')
016.511  000.061  000.061: require('nvim-tree.actions.fs.remove-file')
016.572  000.059  000.059: require('nvim-tree.actions.fs.rename-file')
016.607  000.033  000.033: require('nvim-tree.actions.fs.trash')
016.610  000.367  000.083: require('nvim-tree.actions.fs')
016.683  000.030  000.030: require('nvim-tree.actions.moves.item')
016.714  000.030  000.030: require('nvim-tree.actions.moves.parent')
016.740  000.025  000.025: require('nvim-tree.actions.moves.sibling')
016.742  000.131  000.047: require('nvim-tree.actions.moves')
016.804  000.030  000.030: require('nvim-tree.actions.node.file-popup')
016.875  000.069  000.069: require('nvim-tree.actions.node.open-file')
016.902  000.026  000.026: require('nvim-tree.actions.node.run-command')
016.946  000.043  000.043: require('nvim-tree.actions.node.system-open')
016.948  000.205  000.037: require('nvim-tree.actions.node')
017.037  000.048  000.048: require('nvim-tree.actions.root.change-dir')
017.073  000.035  000.035: require('nvim-tree.actions.root.dir-up')
017.074  000.125  000.041: require('nvim-tree.actions.root')
017.145  000.032  000.032: require('nvim-tree.actions.tree.find-file')
017.216  000.036  000.036: require('nvim-tree.actions.tree.modifiers.collapse-all')
017.257  000.041  000.041: require('nvim-tree.actions.tree.modifiers.expand-all')
017.308  000.049  000.049: require('nvim-tree.actions.tree.modifiers.toggles')
017.309  000.164  000.037: require('nvim-tree.actions.tree.modifiers')
017.349  000.039  000.039: require('nvim-tree.actions.tree.open')
017.393  000.043  000.043: require('nvim-tree.actions.tree.toggle')
017.394  000.319  000.042: require('nvim-tree.actions.tree')
017.395  001.334  000.043: require('nvim-tree.actions')
017.441  000.045  000.045: require('nvim-tree.appearance.diagnostics')
017.560  000.064  000.064: require('nvim-tree.keymap')
017.568  000.125  000.061: require('nvim-tree.help')
017.941  000.372  000.372: require('nvim-tree.marks.navigation')
017.998  000.055  000.055: require('nvim-tree.marks.bulk-delete')
018.030  000.031  000.031: require('nvim-tree.marks.bulk-trash')
018.063  000.032  000.032: require('nvim-tree.marks.bulk-move')
018.096  002.079  000.085: require('nvim-tree.api')
018.110  002.129  000.050: require('nvim-tree.commands')
018.154  000.043  000.043: require('nvim-tree.legacy')
018.192  005.049  000.690: require('nvim-tree')
020.134  000.970  000.970: require('nvim-web-devicons.icons-default')
020.308  001.406  000.436: require('nvim-web-devicons')
022.586  000.189  000.189: sourcing /home/gary/.local/share/nvim/lazy/leap.nvim/plugin/init.lua
023.094  000.024  000.024: require('bufferline.lazy')
023.148  000.051  000.051: require('bufferline.commands')
023.201  000.050  000.050: require('bufferline.config')
023.203  000.374  000.249: require('bufferline')
023.293  000.037  000.037: require('bufferline.utils')
023.294  000.088  000.050: require('bufferline.groups')
023.327  000.024  000.024: require('bufferline.constants')
023.358  000.029  000.029: require('bufferline.colors')
023.444  000.034  000.034: require('bufferline.highlights')
023.946  000.090  000.090: require('vim.version')
024.642  000.821  000.731: require('bufferline.hover')
024.811  000.091  000.091: require('bufferline.ui')
024.918  000.021  000.021: sourcing /home/gary/.local/share/nvim/lazy/vim-numbertoggle/plugin/number_toggle.vim
025.391  000.031  000.031: require('toggleterm.lazy')
025.421  000.028  000.028: require('toggleterm.constants')
025.515  000.093  000.093: require('toggleterm.terminal')
025.522  000.337  000.185: require('toggleterm')
025.622  000.038  000.038: require('toggleterm.colors')
025.668  000.044  000.044: require('toggleterm.utils')
025.674  000.150  000.067: require('toggleterm.config')
026.601  000.083  000.083: require('toggleterm.commandline')
028.463  000.165  000.165: require('mason-core.functional')
028.523  000.058  000.058: require('mason-core.path')
028.588  000.063  000.063: require('mason.settings')
028.609  000.396  000.110: require('mason-core.log')
028.613  000.477  000.082: require('mason-core.EventEmitter')
028.690  000.076  000.076: require('mason-core.optional')
028.864  000.107  000.107: require('mason-core.async')
028.940  000.073  000.073: require('mason-core.async.uv')
028.950  000.259  000.080: require('mason-core.fs')
029.067  000.115  000.115: require('mason-registry.sources')
029.173  000.050  000.050: require('mason-core.functional.data')
029.239  000.064  000.064: require('mason-core.functional.function')
029.255  000.186  000.071: require('mason-core.functional.list')
029.270  001.795  000.682: require('mason-registry')
029.275  002.286  000.491: require('mason-tool-installer')
029.284  002.344  000.059: sourcing /home/gary/.local/share/nvim/lazy/mason-tool-installer.nvim/plugin/mason-tool-installer.lua
029.702  000.097  000.097: require('mason-core.functional.relation')
029.814  000.101  000.101: require('mason-core.functional.logic')
029.848  000.431  000.233: require('mason-core.platform')
029.850  000.524  000.094: require('mason')
030.100  000.114  000.114: require('mason-lspconfig.settings')
030.103  000.252  000.138: require('mason-lspconfig')
030.571  000.077  000.077: require('mason-core.functional.string')
030.589  000.265  000.188: require('mason.api.command')
030.801  000.138  000.138: require('mason-lspconfig.notify')
030.803  000.195  000.057: require('mason-lspconfig.lspconfig_hook')
032.534  000.020  000.020: sourcing /home/gary/.local/share/nvim/lazy/plenary.nvim/plugin/plenary.vim
033.063  000.030  000.030: require('plenary.bit')
033.091  000.026  000.026: require('plenary.functional')
033.116  000.383  000.327: require('plenary.path')
033.162  001.012  000.609: require('plenary.log')
033.194  001.158  000.146: require('lsp-file-operations.log')
033.196  001.536  000.377: require('lsp-file-operations')
034.195  000.155  000.155: require('neodev')
034.224  000.027  000.027: require('neodev.config')
034.302  000.026  000.026: require('neodev.util')
034.304  000.075  000.049: require('neodev.lsp')
034.449  000.038  000.038: require('vim.lsp.log')
034.735  000.284  000.284: require('vim.lsp.protocol')
034.895  000.081  000.081: require('vim.lsp._snippet_grammar')
034.921  000.025  000.025: require('vim.highlight')
034.930  000.192  000.085: require('vim.lsp.util')
034.979  000.022  000.022: require('vim.lsp.sync')
034.982  000.051  000.029: require('vim.lsp._changetracking')
035.026  000.044  000.044: require('vim.lsp.rpc')
035.062  000.723  000.115: require('vim.lsp')
035.109  000.804  000.081: require('lspconfig.util')
035.356  000.069  000.069: sourcing /home/gary/.local/share/nvim/lazy/nvim-lspconfig/plugin/lspconfig.lua
035.486  000.021  000.021: require('lspconfig.async')
035.489  000.053  000.032: require('lspconfig.configs')
035.490  000.112  000.059: require('lspconfig')
035.536  000.022  000.022: require('cmp_nvim_lsp.source')
035.538  000.047  000.025: require('cmp_nvim_lsp')
035.712  000.025  000.025: require('mason-core.functional.table')
035.741  000.184  000.159: require('mason-lspconfig.mappings.server')
035.850  000.030  000.030: require('lspconfig.server_configurations.vimls')
036.011  000.033  000.033: require('lspconfig.manager')
036.040  000.025  000.025: require('lspconfig.server_configurations.lua_ls')
036.183  000.035  000.035: require('lspconfig.server_configurations.clangd')
036.403  000.041  000.041: require('lspconfig.server_configurations.html')
036.497  000.031  000.031: require('lspconfig.server_configurations.pyright')
036.594  000.033  000.033: require('lspconfig.server_configurations.ltex')
036.773  000.033  000.033: require('lspconfig.server_configurations.cssls')
036.925  000.035  000.035: require('lspconfig.server_configurations.bashls')
037.046  006.242  002.937: require('lspconfig.util')
037.079  000.031  000.031: require('mason-lspconfig.server_config_extensions')
037.145  000.041  000.041: require('lspconfig.server_configurations.omnisharp')
037.193  000.030  000.030: require('mason-lspconfig.ensure_installed')
037.270  000.034  000.034: require('mason-core.result')
037.400  000.062  000.062: require('mason-core.process')
037.524  000.122  000.122: require('mason-core.spawn')
037.526  000.224  000.040: require('mason-core.managers.powershell')
037.554  000.027  000.027: require('mason.version')
037.556  000.285  000.034: require('mason-core.fetch')
037.587  000.030  000.030: require('mason-core.providers')
037.755  000.073  000.073: require('mason-core.purl')
037.765  000.135  000.063: require('mason-core.package')
037.916  000.039  000.039: require('mason-core.installer.registry.expr')
037.922  000.091  000.052: require('mason-core.installer.registry.link')
038.164  000.067  000.067: require('mason-core.receipt')
038.180  000.133  000.066: require('mason-core.installer.context')
038.228  000.047  000.047: require('mason-core.async.control')
038.280  000.051  000.051: require('mason-core.installer.linker')
038.284  000.282  000.051: require('mason-core.installer')
038.303  000.350  000.068: require('mason-core.installer.managers.std')
038.305  000.381  000.031: require('mason-core.installer.registry.schemas')
038.349  000.044  000.044: require('mason-core.installer.registry.util')
038.361  000.595  000.079: require('mason-core.installer.registry')
038.363  000.775  000.044: require('mason-registry.sources.util')
038.374  001.175  000.052: require('mason-registry.sources.github')
041.284  000.046  000.046: require('mason-core.functional.number')
041.301  000.118  000.073: require('mason-lspconfig.api.command')
042.041  000.014  000.014: sourcing /home/gary/.local/share/nvim/lazy/todo-comments.nvim/plugin/todo.vim
042.471  000.101  000.101: require('todo-comments.util')
042.481  000.192  000.090: require('todo-comments.config')
042.606  000.082  000.082: require('todo-comments.highlight')
042.612  000.130  000.048: require('todo-comments.jump')
042.614  000.535  000.213: require('todo-comments')
042.919  000.187  000.187: sourcing /home/gary/.local/share/nvim/lazy/telescope.nvim/plugin/telescope.lua
043.075  000.066  000.066: require('telescope._extensions')
043.081  000.140  000.074: require('telescope')
043.524  000.104  000.104: require('plenary.strings')
043.601  000.075  000.075: require('telescope.deprecated')
043.763  000.084  000.084: require('telescope.log')
044.011  000.058  000.058: require('plenary.compat')
044.018  000.162  000.104: require('plenary.job')
044.065  000.046  000.046: require('telescope.state')
044.075  000.310  000.103: require('telescope.utils')
044.081  000.478  000.084: require('telescope.sorters')
045.388  002.137  001.480: require('telescope.config')
045.829  000.312  000.312: require('plenary.window.border')
045.964  000.133  000.133: require('plenary.window')
046.039  000.074  000.074: require('plenary.popup.utils')
046.056  000.666  000.147: require('plenary.popup')
046.153  000.097  000.097: require('telescope.pickers.scroller')
046.244  000.089  000.089: require('telescope.actions.state')
046.335  000.090  000.090: require('telescope.actions.utils')
046.467  000.035  000.035: require('telescope.actions.mt')
046.474  000.137  000.102: require('telescope.actions.set')
046.558  000.037  000.037: require('telescope.config.resolve')
046.560  000.084  000.048: require('telescope.pickers.entry_display')
046.592  000.032  000.032: require('telescope.from_entry')
046.745  003.663  000.331: require('telescope.actions')
047.544  000.158  000.158: require('fzf_lib')
047.547  000.211  000.053: require('telescope._extensions.fzf')
047.974  000.025  000.025: require('auto-session.logger')
047.981  000.066  000.041: require('auto-session.lib')
048.008  000.026  000.026: require('auto-session.autocmds')
048.025  000.325  000.233: require('auto-session')
048.093  000.025  000.025: require('auto-session.session-lens.library')
048.118  000.023  000.023: require('auto-session.session-lens.actions')
048.122  000.088  000.039: require('auto-session.session-lens')
048.722  000.181  000.181: require('leap-spooky')
049.626  000.031  000.031: require('lualine_require')
049.786  000.420  000.388: require('lualine')
049.813  000.026  000.026: require('lazy.status')
050.989  000.030  000.030: require('lualine.utils.mode')
053.992  000.043  000.043: sourcing /usr/share/nvim/runtime/plugin/editorconfig.lua
054.120  000.107  000.107: sourcing /usr/share/nvim/runtime/plugin/gzip.vim
054.226  000.088  000.088: sourcing /usr/share/nvim/runtime/plugin/man.lua
054.607  000.118  000.118: sourcing /usr/share/nvim/runtime/pack/dist/opt/matchit/plugin/matchit.vim
054.617  000.371  000.254: sourcing /usr/share/nvim/runtime/plugin/matchit.vim
054.705  000.069  000.069: sourcing /usr/share/nvim/runtime/plugin/matchparen.vim
054.728  000.007  000.007: sourcing /usr/share/nvim/runtime/plugin/netrwPlugin.vim
054.814  000.067  000.067: sourcing /usr/share/nvim/runtime/plugin/osc52.lua
054.919  000.086  000.086: sourcing /usr/share/nvim/runtime/plugin/rplugin.vim
054.972  000.035  000.035: sourcing /usr/share/nvim/runtime/plugin/shada.vim
055.006  000.010  000.010: sourcing /usr/share/nvim/runtime/plugin/spellfile.vim
055.066  000.042  000.042: sourcing /usr/share/nvim/runtime/plugin/tarPlugin.vim
055.153  000.067  000.067: sourcing /usr/share/nvim/runtime/plugin/tohtml.lua
055.189  000.013  000.013: sourcing /usr/share/nvim/runtime/plugin/tutor.vim
055.291  000.083  000.083: sourcing /usr/share/nvim/runtime/plugin/zipPlugin.vim
055.687  000.118  000.118: require('bigfile.features')
055.690  000.159  000.041: require('bigfile')
055.740  000.250  000.091: sourcing /home/gary/.local/share/nvim/lazy/bigfile.nvim/after/plugin/bigfile.lua
055.819  000.030  000.030: sourcing /home/gary/.local/share/nvim/lazy/cmp-nvim-lsp/after/plugin/cmp_nvim_lsp.lua
055.852  052.251  018.680: require('KaGaSi.lazy')
055.853  053.787  000.010: sourcing /home/gary/.config/nvim/init.lua
055.858  000.277: sourcing vimrc file(s)
056.168  000.047  000.047: sourcing /usr/share/nvim/runtime/filetype.lua
056.346  000.048  000.048: sourcing /usr/share/nvim/runtime/syntax/synload.vim
057.000  000.787  000.739: sourcing /usr/share/nvim/runtime/syntax/syntax.vim
057.008  000.316: inits 3
058.321  001.313: reading ShaDa
058.343  000.023: making windows
059.078  000.072  000.072: require('gitsigns.async')
059.135  000.055  000.055: require('gitsigns.debug.log')
059.303  000.167  000.167: require('gitsigns.config')
059.306  000.624  000.331: require('gitsigns')
059.478  000.112  000.112: require('gitsigns.highlight')
062.844  000.198  000.198: require('gitsigns.util')
062.972  000.122  000.122: require('gitsigns.system')
063.049  000.073  000.073: require('gitsigns.message')
063.121  000.069  000.069: require('gitsigns.git.version')
063.135  003.240  002.777: require('gitsigns.git')
063.346  000.072  000.072: require('gitsigns.cache')
063.389  000.041  000.041: require('gitsigns.signs')
063.429  000.039  000.039: require('gitsigns.status')
063.464  000.034  000.034: require('gitsigns.debounce')
063.499  000.034  000.034: require('gitsigns.diff')
063.855  000.353  000.353: require('gitsigns.hunks')
063.864  000.727  000.154: require('gitsigns.manager')
063.872  004.077  000.110: require('gitsigns.attach')
063.957  000.060  000.060: require('gitsigns.current_line_blame')
064.143  000.056  000.056: require('vim._system')
066.215  000.034  000.034: require('nvim-surround.input')
066.244  000.026  000.026: require('nvim-surround.functional')
066.257  000.172  000.112: require('nvim-surround.config')
066.264  000.234  000.062: require('nvim-surround.buffer')
066.293  000.028  000.028: require('nvim-surround.cache')
066.349  000.055  000.055: require('nvim-surround.utils')
066.352  000.685  000.367: require('nvim-surround')
067.225  000.059  000.059: require('ibl.utils')
067.233  000.113  000.054: require('ibl.config')
067.297  000.027  000.027: require('ibl.indent')
067.302  000.068  000.041: require('ibl.hooks')
067.303  000.216  000.035: require('ibl.highlights')
067.528  000.224  000.224: require('ibl.autocmds')
067.866  000.337  000.337: require('ibl.inlay_hints')
067.956  000.066  000.066: require('ibl.virt_text')
068.376  000.360  000.360: require('ibl.scope_languages')
068.378  000.421  000.061: require('ibl.scope')
068.386  001.546  000.282: require('ibl')
068.404  001.601  000.055: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/plugin/commands.lua
068.680  000.235  000.235: require('vim.iter')
068.958  000.261  000.261: require('vim.lsp.handlers')
069.612  000.138  000.138: require('vim.treesitter.language')
069.664  000.050  000.050: require('vim.func')
069.754  000.087  000.087: require('vim.func._memoize')
069.774  000.479  000.204: require('vim.treesitter.query')
069.842  000.067  000.067: require('vim.treesitter._range')
069.856  000.669  000.124: require('vim.treesitter.languagetree')
069.867  000.779  000.110: require('vim.treesitter')
071.341  000.149  000.149: require('nvim-treesitter.utils')
072.854  001.511  001.511: require('nvim-treesitter.parsers')
073.193  000.069  000.069: require('nvim-treesitter.compat')
073.386  000.132  000.132: require('nvim-treesitter.ts_utils')
073.392  000.197  000.066: require('nvim-treesitter.tsrange')
073.439  000.046  000.046: require('nvim-treesitter.caching')
073.451  000.432  000.120: require('nvim-treesitter.query')
073.469  000.549  000.117: require('nvim-treesitter.configs')
073.472  000.616  000.067: require('nvim-treesitter.info')
073.574  000.101  000.101: require('nvim-treesitter.shell_command_selectors')
073.638  002.582  000.205: require('nvim-treesitter.install')
073.692  000.053  000.053: require('nvim-treesitter.statusline')
073.760  000.068  000.068: require('nvim-treesitter.query_predicates')
073.762  003.178  000.475: require('nvim-treesitter')
073.943  003.409  000.231: sourcing /home/gary/.local/share/nvim/lazy/nvim-treesitter/plugin/nvim-treesitter.lua
075.838  000.214  000.214: require('nvim-treesitter.indent')
076.252  000.113  000.113: require('nvim-treesitter.highlight')
076.704  000.211  000.211: require('nvim-treesitter.locals')
076.713  000.320  000.109: require('nvim-treesitter.incremental_selection')
077.200  000.161  000.161: sourcing /home/gary/.local/share/nvim/lazy/nvim-ts-context-commentstring/plugin/ts_context_commentstring.lua
078.195  000.232  000.232: require('Comment.config')
078.985  000.592  000.592: require('Comment.ft')
079.002  000.805  000.213: require('Comment.utils')
079.150  000.147  000.147: require('Comment.opfunc')
079.271  000.121  000.121: require('Comment.extra')
079.287  001.978  000.674: require('Comment.api')
079.362  002.102  000.124: sourcing /home/gary/.local/share/nvim/lazy/Comment.nvim/plugin/Comment.lua
079.413  000.034  000.034: require('Comment')
079.446  000.032  000.032: require('ts_context_commentstring.integrations.comment_nvim')
082.916  000.190  000.190: sourcing /usr/share/nvim/runtime/autoload/provider/clipboard.vim
084.600  000.471  000.471: require('vim.filetype.detect')
085.221  000.057  000.057: sourcing /usr/share/nvim/runtime/ftplugin/c.vim
085.322  000.088  000.088: sourcing /usr/share/nvim/runtime/ftplugin/c.lua
085.672  000.114  000.114: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/ftplugin/c.lua
086.328  000.026  000.026: sourcing /usr/share/nvim/runtime/indent/c.vim
087.304  000.153  000.153: require('vim.lsp.client')
087.457  000.042  000.042: require('vim.glob')
087.461  000.091  000.049: require('vim.lsp._dynamic')
092.404  000.642  000.642: sourcing /usr/share/nvim/runtime/syntax/c.vim
128.403  000.775  000.775: require('vim.treesitter.highlighter')
129.521  000.058  000.058: require('ts_context_commentstring.utils')
129.601  000.078  000.078: require('ts_context_commentstring.config')
129.606  000.212  000.076: require('ts_context_commentstring.internal')
129.743  000.074  000.074: require('editorconfig')
129.863  053.749: opening buffers
129.929  000.045  000.045: require('bufferline.state')
131.963  002.055: BufEnter autocommands
136.108  000.037  000.037: sourcing /usr/share/nvim/runtime/ftplugin/conf.vim
137.374  000.037  000.037: sourcing /usr/share/nvim/runtime/syntax/conf.vim
146.303  014.267: editing files in windows
151.551  005.248: executing command arguments
152.332  000.558  000.558: require('alpha')
152.399  000.064  000.064: require('alpha.themes.dashboard')
166.060  013.887: VimEnter autocommands
166.097  000.036: UIEnter autocommands
166.105  000.008: diff scrollbinding
166.106  000.001: before starting main loop
166.557  000.121  000.121: require('bufferline.tabpages')
166.713  000.113  000.113: require('bufferline.models')
166.774  000.059  000.059: require('bufferline.pick')
166.828  000.052  000.052: require('bufferline.duplicates')
166.951  000.116  000.116: require('bufferline.diagnostics')
167.281  000.077  000.077: require('bufferline.numbers')
167.409  000.080  000.080: require('bufferline.sorters')
167.523  000.076  000.076: require('bufferline.offset')
167.602  000.076  000.076: require('bufferline.custom_area')
208.100  041.226: first screen update
208.106  000.006: --- NVIM STARTED ---

--- Startup times for process: Primary/TUI ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.164  000.162: event init
000.230  000.066: early init
000.267  000.037: locale set
000.300  000.033: init first window
000.495  000.196: inits 1
000.502  000.006: window checked
000.565  000.064: parsing arguments
000.929  000.029  000.029: require('vim.shared')
000.989  000.030  000.030: require('vim.inspect')
001.029  000.031  000.031: require('vim._options')
001.030  000.099  000.038: require('vim._editor')
001.031  000.144  000.016: require('vim._init_packages')
001.032  000.323: init lua interpreter
001.652  000.619: --- NVIM STARTED ---

--- Startup times for process: Embedded ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.123  000.122: event init
000.182  000.059: early init
000.241  000.058: locale set
000.314  000.073: init first window
000.512  000.198: inits 1
000.522  000.011: window checked
000.580  000.058: parsing arguments
000.999  000.033  000.033: require('vim.shared')
001.231  000.156  000.156: require('vim.inspect')
001.306  000.058  000.058: require('vim._options')
001.307  000.306  000.093: require('vim._editor')
001.308  000.375  000.036: require('vim._init_packages')
001.310  000.355: init lua interpreter
001.536  000.226: expanding arguments
001.580  000.043: inits 2
001.881  000.302: init highlight
001.882  000.001: waiting for UI
002.210  000.328: done waiting for UI
002.213  000.003: clear screen
002.289  000.005  000.005: require('vim.keymap')
002.694  000.478  000.474: require('vim._defaults')
002.696  000.004: init default mappings & autocommands
003.046  000.096  000.096: sourcing /usr/share/nvim/runtime/ftplugin.vim
003.104  000.022  000.022: sourcing /usr/share/nvim/runtime/indent.vim
003.250  000.095  000.095: sourcing /usr/share/nvim/archlinux.lua
003.254  000.126  000.032: sourcing /etc/xdg/nvim/sysinit.vim
004.165  000.645  000.645: require('KaGaSi.core.keymap')
004.980  000.813  000.813: require('KaGaSi.core.options')
005.053  000.071  000.071: require('KaGaSi.core.fold')
005.085  000.030  000.030: require('KaGaSi.core.filepos')
005.105  000.019  000.019: require('KaGaSi.core.pdf')
005.106  001.809  000.231: require('KaGaSi.core')
005.490  000.289  000.289: require('lazy')
005.515  000.012  000.012: require('ffi')
005.552  000.018  000.018: require('vim.fs')
005.712  000.126  000.126: require('vim.uri')
005.720  000.202  000.059: require('vim.loader')
005.867  000.135  000.135: require('lazy.stats')
005.951  000.067  000.067: require('lazy.core.util')
006.017  000.064  000.064: require('lazy.core.config')
006.191  000.075  000.075: require('lazy.core.handler')
006.385  000.193  000.193: require('lazy.core.plugin')
006.410  000.392  000.124: require('lazy.core.loader')
007.351  000.169  000.169: require('KaGaSi.plugins.alpha')
007.518  000.093  000.093: require('KaGaSi.plugins.auto-session')
007.590  000.053  000.053: require('KaGaSi.plugins.autopairs')
007.727  000.117  000.117: require('KaGaSi.plugins.bigfile')
007.908  000.169  000.169: require('KaGaSi.plugins.bufferline')
007.987  000.042  000.042: require('KaGaSi.plugins.colorscheme')
008.034  000.032  000.032: require('KaGaSi.plugins.comment')
008.095  000.034  000.034: require('KaGaSi.plugins.dressing')
008.132  000.029  000.029: require('KaGaSi.plugins.gitsigns')
008.182  000.045  000.045: require('KaGaSi.plugins.indent-blankline')
008.223  000.036  000.036: require('KaGaSi.plugins.leap')
008.257  000.029  000.029: require('KaGaSi.plugins.leap-spooky')
008.316  000.048  000.048: require('KaGaSi.plugins.lspconfig')
008.602  000.100  000.100: require('KaGaSi.plugins.ltex_extra')
008.666  000.060  000.060: require('KaGaSi.plugins.lualine')
008.712  000.024  000.024: require('KaGaSi.plugins.mason')
008.755  000.022  000.022: require('KaGaSi.plugins.numbertoggle')
008.783  000.023  000.023: require('KaGaSi.plugins.nvim-cmp')
008.860  000.057  000.057: require('KaGaSi.plugins.nvim-tree')
008.896  000.026  000.026: require('KaGaSi.plugins.surround')
008.931  000.028  000.028: require('KaGaSi.plugins.telescope')
008.975  000.026  000.026: require('KaGaSi.plugins.todo-comments')
009.021  000.028  000.028: require('KaGaSi.plugins.toggleterm')
009.103  000.074  000.074: require('KaGaSi.plugins.treesitter')
009.138  000.030  000.030: require('KaGaSi.plugins.trouble')
009.546  000.028  000.028: require('KaGaSi.plugins.vim-maximizer')
009.573  000.022  000.022: require('KaGaSi.plugins.vimtex')
010.591  000.950  000.950: require('vim.filetype')
010.594  001.015  000.065: require('KaGaSi.plugins.vimwiki')
010.627  000.028  000.028: require('KaGaSi.plugins.which_key')
010.858  000.033  000.033: require('lazy.core.handler.cmd')
010.889  000.029  000.029: require('lazy.core.handler.event')
010.923  000.030  000.030: require('lazy.core.handler.ft')
010.963  000.039  000.039: require('lazy.core.handler.keys')
011.259  000.016  000.016: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/cls.vim
011.285  000.009  000.009: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tex.vim
011.310  000.008  000.008: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tikz.vim
011.719  000.119  000.119: sourcing /usr/share/nvim/runtime/filetype.lua
011.961  000.107  000.107: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/scroll.vim
012.010  000.028  000.028: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/themer.lua
012.943  000.246  000.246: require('themer')
012.995  000.049  000.049: require('themer.config')
013.072  000.033  000.033: require('themer.modules.core')
013.101  000.027  000.027: require('themer.modules.core.api')
013.152  000.050  000.050: require('themer.modules.themes.sonokai_deep')
013.187  000.033  000.033: require('themer.modules.core.utils')
013.293  000.105  000.105: require('themer.modules.core.mapper')
014.577  001.906  001.363: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/colors/themer_sonokai_deep.lua
014.833  000.024  000.024: sourcing /home/gary/.local/share/nvim/lazy/nvim-web-devicons/plugin/nvim-web-devicons.vim
015.261  000.039  000.039: require('bufferline.lazy')
015.366  000.095  000.095: require('bufferline.commands')
015.512  000.143  000.143: require('bufferline.config')
015.518  000.624  000.347: require('bufferline')
015.713  000.072  000.072: require('bufferline.utils')
015.715  000.193  000.121: require('bufferline.groups')
015.782  000.049  000.049: require('bufferline.constants')
015.828  000.043  000.043: require('bufferline.colors')
016.109  000.072  000.072: require('bufferline.highlights')
016.234  000.005  000.005: require('vim.F')
016.816  000.172  000.172: require('vim.version')
017.923  001.359  001.187: require('bufferline.hover')
018.150  000.134  000.134: require('bufferline.ui')
018.551  000.043  000.043: require('toggleterm.lazy')
018.594  000.041  000.041: require('toggleterm.constants')
018.692  000.096  000.096: require('toggleterm.terminal')
018.696  000.355  000.175: require('toggleterm')
018.789  000.036  000.036: require('toggleterm.colors')
018.826  000.036  000.036: require('toggleterm.utils')
018.829  000.132  000.061: require('toggleterm.config')
019.530  000.064  000.064: require('toggleterm.commandline')
020.008  000.229  000.229: sourcing /home/gary/.local/share/nvim/lazy/leap.nvim/plugin/init.lua
020.434  000.318  000.318: require('leap-spooky')
021.522  000.029  000.029: require('auto-session.logger')
021.530  000.080  000.051: require('auto-session.lib')
021.560  000.029  000.029: require('auto-session.autocmds')
021.578  000.415  000.306: require('auto-session')
021.866  000.015  000.015: sourcing /home/gary/.local/share/nvim/lazy/plenary.nvim/plugin/plenary.vim
022.113  000.011  000.011: sourcing /home/gary/.local/share/nvim/lazy/todo-comments.nvim/plugin/todo.vim
022.424  000.036  000.036: require('todo-comments.util')
022.433  000.089  000.053: require('todo-comments.config')
022.565  000.081  000.081: require('todo-comments.highlight')
022.575  000.141  000.060: require('todo-comments.jump')
022.577  000.454  000.224: require('todo-comments')
022.790  000.163  000.163: sourcing /home/gary/.local/share/nvim/lazy/telescope.nvim/plugin/telescope.lua
022.990  000.075  000.075: require('telescope._extensions')
022.994  000.190  000.115: require('telescope')
023.659  000.067  000.067: require('plenary.bit')
023.749  000.087  000.087: require('plenary.functional')
023.806  000.301  000.147: require('plenary.path')
023.818  000.395  000.094: require('plenary.strings')
023.934  000.114  000.114: require('telescope.deprecated')
024.641  000.343  000.343: require('plenary.log')
024.702  000.637  000.293: require('telescope.log')
024.968  000.034  000.034: require('plenary.compat')
024.974  000.123  000.089: require('plenary.job')
025.100  000.093  000.093: require('telescope.state')
025.114  000.410  000.193: require('telescope.utils')
025.120  001.179  000.133: require('telescope.sorters')
026.266  002.997  001.309: require('telescope.config')
026.356  000.035  000.035: require('plenary.window.border')
026.383  000.026  000.026: require('plenary.window')
026.405  000.021  000.021: require('plenary.popup.utils')
026.408  000.139  000.057: require('plenary.popup')
026.489  000.081  000.081: require('telescope.pickers.scroller')
026.553  000.062  000.062: require('telescope.actions.state')
026.595  000.041  000.041: require('telescope.actions.utils')
026.667  000.034  000.034: require('telescope.actions.mt')
026.679  000.083  000.048: require('telescope.actions.set')
026.750  000.036  000.036: require('telescope.config.resolve')
026.752  000.073  000.037: require('telescope.pickers.entry_display')
026.783  000.031  000.031: require('telescope.from_entry')
026.857  003.862  000.357: require('telescope.actions')
027.554  000.098  000.098: require('fzf_lib')
027.557  000.146  000.048: require('telescope._extensions.fzf')
027.601  005.985  001.143: require('telescope')
027.666  000.024  000.024: require('auto-session.session-lens.library')
027.691  000.025  000.025: require('auto-session.session-lens.actions')
027.719  000.117  000.068: require('auto-session.session-lens')
028.542  000.041  000.041: require('mason-core.functional')
028.567  000.023  000.023: require('mason-core.path')
028.593  000.025  000.025: require('mason.settings')
028.603  000.126  000.037: require('mason-core.log')
028.605  000.149  000.022: require('mason-core.EventEmitter')
028.628  000.023  000.023: require('mason-core.optional')
028.696  000.036  000.036: require('mason-core.async')
028.726  000.029  000.029: require('mason-core.async.uv')
028.731  000.101  000.037: require('mason-core.fs')
028.766  000.035  000.035: require('mason-registry.sources')
028.837  000.027  000.027: require('mason-core.functional.data')
028.868  000.029  000.029: require('mason-core.functional.function')
028.878  000.110  000.054: require('mason-core.functional.list')
028.887  000.493  000.075: require('mason-registry')
028.892  000.674  000.181: require('mason-tool-installer')
028.900  000.709  000.035: sourcing /home/gary/.local/share/nvim/lazy/mason-tool-installer.nvim/plugin/mason-tool-installer.lua
029.042  000.028  000.028: require('mason-core.functional.relation')
029.076  000.029  000.029: require('mason-core.functional.logic')
029.084  000.121  000.064: require('mason-core.platform')
029.085  000.156  000.035: require('mason')
029.155  000.034  000.034: require('mason-lspconfig.settings')
029.157  000.071  000.037: require('mason-lspconfig')
029.339  000.034  000.034: require('mason-core.functional.string')
029.353  000.102  000.068: require('mason.api.command')
029.446  000.024  000.024: require('mason-lspconfig.notify')
029.449  000.062  000.037: require('mason-lspconfig.lspconfig_hook')
030.292  000.076  000.076: require('lsp-file-operations.log')
030.294  000.286  000.209: require('lsp-file-operations')
031.065  000.044  000.044: require('nvim-tree.notify')
031.071  000.095  000.051: require('nvim-tree.events')
031.219  000.055  000.055: require('nvim-tree.log')
031.335  000.053  000.053: require('nvim-tree.iterators.node-iterator')
031.360  000.140  000.087: require('nvim-tree.utils')
031.407  000.046  000.046: require('nvim-tree.git.utils')
031.470  000.061  000.061: require('nvim-tree.git.runner')
031.540  000.069  000.069: require('nvim-tree.watcher')
031.585  000.044  000.044: require('nvim-tree.explorer.node')
031.589  000.479  000.063: require('nvim-tree.git')
031.633  000.043  000.043: require('nvim-tree.explorer.watch')
031.719  000.047  000.047: require('nvim-tree.explorer.node-builders')
031.779  000.059  000.059: require('nvim-tree.explorer.sorters')
031.871  000.044  000.044: require('nvim-tree.marks')
031.876  000.096  000.052: require('nvim-tree.explorer.filters')
032.075  000.150  000.150: require('nvim-tree.view')
032.081  000.205  000.054: require('nvim-tree.live-filter')
032.083  000.449  000.043: require('nvim-tree.explorer.explore')
032.144  000.060  000.060: require('nvim-tree.explorer.reload')
032.147  001.075  000.044: require('nvim-tree.explorer')
032.149  001.204  000.033: require('nvim-tree.core')
032.221  000.071  000.071: require('nvim-tree.renderer.components.padding')
032.280  000.057  000.057: require('nvim-tree.renderer.components.icons')
032.340  000.059  000.059: require('nvim-tree.renderer.components.full-name')
032.506  000.045  000.045: require('nvim-tree.enum')
032.550  000.043  000.043: require('nvim-tree.renderer.decorator')
032.552  000.144  000.056: require('nvim-tree.renderer.decorator.bookmarks')
032.593  000.039  000.039: require('nvim-tree.renderer.decorator.copied')
032.633  000.039  000.039: require('nvim-tree.renderer.decorator.cut')
032.750  000.061  000.061: require('nvim-tree.diagnostics')
033.226  000.474  000.474: require('vim.diagnostic')
033.232  000.599  000.063: require('nvim-tree.renderer.decorator.diagnostics')
033.386  000.153  000.153: require('nvim-tree.renderer.decorator.git')
033.665  000.136  000.136: require('nvim-tree.buffers')
033.667  000.280  000.144: require('nvim-tree.renderer.decorator.modified')
033.779  000.111  000.111: require('nvim-tree.renderer.decorator.opened')
033.794  001.453  000.087: require('nvim-tree.renderer.builder')
033.801  002.935  000.091: require('nvim-tree.renderer')
033.812  003.045  000.110: require('nvim-tree.lib')
033.985  000.173  000.173: require('nvim-tree.appearance')
034.149  000.024  000.024: require('nvim-tree.actions.finders.find-file')
034.198  000.048  000.048: require('nvim-tree.actions.finders.search-node')
034.199  000.097  000.025: require('nvim-tree.actions.finders')
034.279  000.024  000.024: require('nvim-tree.actions.reloaders')
034.290  000.066  000.042: require('nvim-tree.actions.fs.copy-paste')
034.315  000.024  000.024: require('nvim-tree.actions.fs.create-file')
034.354  000.038  000.038: require('nvim-tree.actions.fs.remove-file')
034.399  000.043  000.043: require('nvim-tree.actions.fs.rename-file')
034.439  000.038  000.038: require('nvim-tree.actions.fs.trash')
034.441  000.241  000.031: require('nvim-tree.actions.fs')
034.539  000.033  000.033: require('nvim-tree.actions.moves.item')
034.563  000.022  000.022: require('nvim-tree.actions.moves.parent')
034.585  000.021  000.021: require('nvim-tree.actions.moves.sibling')
034.586  000.145  000.068: require('nvim-tree.actions.moves')
034.674  000.036  000.036: require('nvim-tree.actions.node.file-popup')
034.747  000.071  000.071: require('nvim-tree.actions.node.open-file')
034.775  000.027  000.027: require('nvim-tree.actions.node.run-command')
034.799  000.023  000.023: require('nvim-tree.actions.node.system-open')
034.800  000.213  000.056: require('nvim-tree.actions.node')
034.857  000.033  000.033: require('nvim-tree.actions.root.change-dir')
034.888  000.030  000.030: require('nvim-tree.actions.root.dir-up')
034.890  000.089  000.025: require('nvim-tree.actions.root')
034.952  000.029  000.029: require('nvim-tree.actions.tree.find-file')
035.025  000.030  000.030: require('nvim-tree.actions.tree.modifiers.collapse-all')
035.057  000.031  000.031: require('nvim-tree.actions.tree.modifiers.expand-all')
035.096  000.038  000.038: require('nvim-tree.actions.tree.modifiers.toggles')
035.097  000.144  000.045: require('nvim-tree.actions.tree.modifiers')
035.129  000.031  000.031: require('nvim-tree.actions.tree.open')
035.161  000.031  000.031: require('nvim-tree.actions.tree.toggle')
035.163  000.273  000.037: require('nvim-tree.actions.tree')
035.164  001.087  000.031: require('nvim-tree.actions')
035.202  000.037  000.037: require('nvim-tree.appearance.diagnostics')
035.259  000.028  000.028: require('nvim-tree.keymap')
035.262  000.059  000.031: require('nvim-tree.help')
035.286  000.023  000.023: require('nvim-tree.marks.navigation')
035.309  000.022  000.022: require('nvim-tree.marks.bulk-delete')
035.331  000.021  000.021: require('nvim-tree.marks.bulk-trash')
035.353  000.021  000.021: require('nvim-tree.marks.bulk-move')
035.446  001.416  000.145: require('nvim-tree.api')
035.452  001.466  000.051: require('nvim-tree.commands')
035.488  000.034  000.034: require('nvim-tree.legacy')
035.538  005.059  000.341: require('nvim-tree')
036.422  000.277  000.277: require('nvim-web-devicons.icons-default')
036.532  000.456  000.179: require('nvim-web-devicons')
038.189  007.891  002.376: require('nvim-tree.api')
038.918  000.168  000.168: require('neodev')
038.945  000.024  000.024: require('neodev.config')
038.995  000.023  000.023: require('neodev.util')
038.997  000.047  000.024: require('neodev.lsp')
039.148  000.058  000.058: require('vim.lsp.log')
039.370  000.220  000.220: require('vim.lsp.protocol')
039.650  000.129  000.129: require('vim.lsp._snippet_grammar')
039.687  000.035  000.035: require('vim.highlight')
039.701  000.328  000.164: require('vim.lsp.util')
039.778  000.036  000.036: require('vim.lsp.sync')
039.782  000.079  000.043: require('vim.lsp._changetracking')
039.856  000.072  000.072: require('vim.lsp.rpc')
039.902  000.872  000.113: require('vim.lsp')
039.959  000.962  000.090: require('lspconfig.util')
040.191  000.058  000.058: sourcing /home/gary/.local/share/nvim/lazy/nvim-lspconfig/plugin/lspconfig.lua
040.318  000.027  000.027: require('lspconfig.async')
040.320  000.073  000.046: require('lspconfig.configs')
040.322  000.109  000.036: require('lspconfig')
040.389  000.032  000.032: require('cmp_nvim_lsp.source')
040.391  000.069  000.037: require('cmp_nvim_lsp')
040.567  000.080  000.080: require('mason-core.functional.table')
040.597  000.186  000.106: require('mason-lspconfig.mappings.server')
040.790  000.058  000.058: require('lspconfig.server_configurations.clangd')
041.086  000.062  000.062: require('lspconfig.manager')
041.163  000.069  000.069: require('lspconfig.server_configurations.pyright')
041.343  000.051  000.051: require('lspconfig.server_configurations.cssls')
041.520  000.051  000.051: require('lspconfig.server_configurations.html')
041.749  000.054  000.054: require('lspconfig.server_configurations.vimls')
041.924  000.068  000.068: require('lspconfig.server_configurations.lua_ls')
042.073  000.071  000.071: require('lspconfig.server_configurations.ltex')
042.316  000.062  000.062: require('lspconfig.server_configurations.bashls')
042.413  012.963  002.616: require('lspconfig.util')
042.466  000.051  000.051: require('mason-lspconfig.server_config_extensions')
042.530  000.062  000.062: require('lspconfig.server_configurations.omnisharp')
042.591  000.044  000.044: require('mason-lspconfig.ensure_installed')
042.731  000.068  000.068: require('mason-core.result')
042.917  000.092  000.092: require('mason-core.process')
043.125  000.206  000.206: require('mason-core.spawn')
043.130  000.353  000.055: require('mason-core.managers.powershell')
043.164  000.032  000.032: require('mason.version')
043.170  000.437  000.052: require('mason-core.fetch')
043.239  000.068  000.068: require('mason-core.providers')
043.485  000.146  000.146: require('mason-core.purl')
043.497  000.221  000.075: require('mason-core.package')
043.710  000.081  000.081: require('mason-core.installer.registry.expr')
043.720  000.164  000.083: require('mason-core.installer.registry.link')
044.469  000.156  000.156: require('mason-core.receipt')
044.494  000.590  000.434: require('mason-core.installer.context')
044.553  000.058  000.058: require('mason-core.async.control')
044.607  000.053  000.053: require('mason-core.installer.linker')
044.612  000.783  000.082: require('mason-core.installer')
044.627  000.862  000.079: require('mason-core.installer.managers.std')
044.628  000.907  000.045: require('mason-core.installer.registry.schemas')
044.672  000.044  000.044: require('mason-core.installer.registry.util')
044.681  001.183  000.068: require('mason-core.installer.registry')
044.683  001.443  000.038: require('mason-registry.sources.util')
044.690  002.092  000.076: require('mason-registry.sources.github')
047.824  000.090  000.090: require('mason-core.functional.number')
047.848  000.233  000.143: require('mason-lspconfig.api.command')
048.663  000.020  000.020: sourcing /home/gary/.local/share/nvim/lazy/vim-numbertoggle/plugin/number_toggle.vim
049.775  000.497  000.497: require('lualine_require')
050.108  001.250  000.753: require('lualine')
050.171  000.062  000.062: require('lazy.status')
051.260  000.030  000.030: require('lualine.utils.mode')
055.679  000.050  000.050: sourcing /usr/share/nvim/runtime/plugin/editorconfig.lua
055.807  000.105  000.105: sourcing /usr/share/nvim/runtime/plugin/gzip.vim
055.889  000.049  000.049: sourcing /usr/share/nvim/runtime/plugin/man.lua
056.385  000.206  000.206: sourcing /usr/share/nvim/runtime/pack/dist/opt/matchit/plugin/matchit.vim
056.432  000.522  000.316: sourcing /usr/share/nvim/runtime/plugin/matchit.vim
056.665  000.176  000.176: sourcing /usr/share/nvim/runtime/plugin/matchparen.vim
056.722  000.016  000.016: sourcing /usr/share/nvim/runtime/plugin/netrwPlugin.vim
056.936  000.175  000.175: sourcing /usr/share/nvim/runtime/plugin/osc52.lua
057.158  000.182  000.182: sourcing /usr/share/nvim/runtime/plugin/rplugin.vim
057.244  000.050  000.050: sourcing /usr/share/nvim/runtime/plugin/shada.vim
057.277  000.009  000.009: sourcing /usr/share/nvim/runtime/plugin/spellfile.vim
057.351  000.053  000.053: sourcing /usr/share/nvim/runtime/plugin/tarPlugin.vim
057.435  000.061  000.061: sourcing /usr/share/nvim/runtime/plugin/tohtml.lua
057.469  000.012  000.012: sourcing /usr/share/nvim/runtime/plugin/tutor.vim
057.567  000.080  000.080: sourcing /usr/share/nvim/runtime/plugin/zipPlugin.vim
057.929  000.123  000.123: require('bigfile.features')
057.933  000.182  000.059: require('bigfile')
057.945  000.243  000.061: sourcing /home/gary/.local/share/nvim/lazy/bigfile.nvim/after/plugin/bigfile.lua
058.046  000.054  000.054: sourcing /home/gary/.local/share/nvim/lazy/cmp-nvim-lsp/after/plugin/cmp_nvim_lsp.lua
058.078  052.971  017.137: require('KaGaSi.lazy')
058.079  054.797  000.017: sourcing /home/gary/.config/nvim/init.lua
058.084  000.347: sourcing vimrc file(s)
058.376  000.029  000.029: sourcing /usr/share/nvim/runtime/filetype.lua
058.542  000.047  000.047: sourcing /usr/share/nvim/runtime/syntax/synload.vim
059.184  000.769  000.722: sourcing /usr/share/nvim/runtime/syntax/syntax.vim
059.193  000.311: inits 3
060.543  001.349: reading ShaDa
060.566  000.023: making windows
061.149  000.042  000.042: require('gitsigns.async')
061.194  000.042  000.042: require('gitsigns.debug.log')
061.308  000.113  000.113: require('gitsigns.config')
061.311  000.454  000.257: require('gitsigns')
061.421  000.078  000.078: require('gitsigns.highlight')
061.839  000.055  000.055: require('gitsigns.util')
061.871  000.031  000.031: require('gitsigns.system')
061.900  000.027  000.027: require('gitsigns.message')
061.928  000.028  000.028: require('gitsigns.git.version')
061.934  000.228  000.087: require('gitsigns.git')
062.031  000.038  000.038: require('gitsigns.cache')
062.064  000.032  000.032: require('gitsigns.signs')
062.090  000.025  000.025: require('gitsigns.status')
062.120  000.028  000.028: require('gitsigns.debounce')
062.144  000.022  000.022: require('gitsigns.diff')
062.192  000.048  000.048: require('gitsigns.hunks')
062.199  000.265  000.072: require('gitsigns.manager')
062.205  000.561  000.068: require('gitsigns.attach')
062.265  000.046  000.046: require('gitsigns.current_line_blame')
062.347  000.042  000.042: require('vim._system')
064.210  000.050  000.050: require('nvim-treesitter.utils')
064.694  000.037  000.037: require('vim.treesitter.language')
064.719  000.023  000.023: require('vim.func')
064.747  000.025  000.025: require('vim.func._memoize')
064.760  000.174  000.088: require('vim.treesitter.query')
064.790  000.029  000.029: require('vim.treesitter._range')
064.801  000.270  000.068: require('vim.treesitter.languagetree')
064.806  000.326  000.056: require('vim.treesitter')
065.095  000.883  000.556: require('nvim-treesitter.parsers')
065.249  000.027  000.027: require('nvim-treesitter.compat')
065.320  000.040  000.040: require('nvim-treesitter.ts_utils')
065.345  000.095  000.054: require('nvim-treesitter.tsrange')
065.376  000.030  000.030: require('nvim-treesitter.caching')
065.382  000.198  000.047: require('nvim-treesitter.query')
065.390  000.257  000.059: require('nvim-treesitter.configs')
065.393  000.297  000.040: require('nvim-treesitter.info')
065.440  000.047  000.047: require('nvim-treesitter.shell_command_selectors')
065.452  001.394  000.117: require('nvim-treesitter.install')
065.480  000.026  000.026: require('nvim-treesitter.statusline')
065.515  000.034  000.034: require('nvim-treesitter.query_predicates')
065.516  001.756  000.302: require('nvim-treesitter')
065.779  002.071  000.314: sourcing /home/gary/.local/share/nvim/lazy/nvim-treesitter/plugin/nvim-treesitter.lua
067.209  000.052  000.052: require('nvim-treesitter.highlight')
067.640  000.065  000.065: require('nvim-treesitter.locals')
067.646  000.134  000.068: require('nvim-treesitter.incremental_selection')
067.796  000.055  000.055: require('nvim-treesitter.indent')
068.572  000.063  000.063: require('ibl.utils')
068.576  000.126  000.063: require('ibl.config')
068.657  000.030  000.030: require('ibl.indent')
068.662  000.085  000.056: require('ibl.hooks')
068.664  000.262  000.051: require('ibl.highlights')
068.701  000.036  000.036: require('ibl.autocmds')
068.739  000.037  000.037: require('ibl.inlay_hints')
068.793  000.053  000.053: require('ibl.virt_text')
069.010  000.185  000.185: require('ibl.scope_languages')
069.012  000.218  000.033: require('ibl.scope')
069.016  000.981  000.375: require('ibl')
069.038  001.039  000.058: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/plugin/commands.lua
069.208  000.083  000.083: require('vim.iter')
069.515  000.244  000.244: require('vim.lsp.handlers')
071.068  000.030  000.030: require('nvim-surround.input')
071.101  000.031  000.031: require('nvim-surround.functional')
071.110  000.171  000.110: require('nvim-surround.config')
071.121  000.232  000.061: require('nvim-surround.buffer')
071.149  000.027  000.027: require('nvim-surround.cache')
071.190  000.040  000.040: require('nvim-surround.utils')
071.193  000.640  000.341: require('nvim-surround')
072.063  000.198  000.198: sourcing /home/gary/.local/share/nvim/lazy/nvim-ts-context-commentstring/plugin/ts_context_commentstring.lua
072.838  000.062  000.062: require('Comment.config')
073.187  000.280  000.280: require('Comment.ft')
073.197  000.357  000.077: require('Comment.utils')
073.265  000.068  000.068: require('Comment.opfunc')
073.365  000.098  000.098: require('Comment.extra')
073.372  001.170  000.585: require('Comment.api')
073.433  001.272  000.102: sourcing /home/gary/.local/share/nvim/lazy/Comment.nvim/plugin/Comment.lua
073.519  000.060  000.060: require('Comment')
073.571  000.050  000.050: require('ts_context_commentstring.integrations.comment_nvim')
076.343  000.189  000.189: sourcing /usr/share/nvim/runtime/autoload/provider/clipboard.vim
078.610  000.946  000.946: require('vim.filetype.detect')
079.221  000.054  000.054: sourcing /usr/share/nvim/runtime/ftplugin/c.vim
079.295  000.060  000.060: sourcing /usr/share/nvim/runtime/ftplugin/c.lua
079.587  000.067  000.067: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/ftplugin/c.lua
079.884  000.012  000.012: sourcing /usr/share/nvim/runtime/indent/c.vim
080.549  000.244  000.244: require('vim.lsp.client')
080.773  000.086  000.086: require('vim.glob')
080.787  000.168  000.082: require('vim.lsp._dynamic')
086.739  000.604  000.604: sourcing /usr/share/nvim/runtime/syntax/c.vim
115.413  000.297  000.297: require('vim.treesitter.highlighter')
128.275  000.032  000.032: require('ts_context_commentstring.utils')
128.316  000.035  000.035: require('ts_context_commentstring.config')
128.318  000.126  000.059: require('ts_context_commentstring.internal')
128.432  000.047  000.047: require('editorconfig')
128.555  058.098: opening buffers
128.618  000.036  000.036: require('bufferline.state')
130.543  001.952: BufEnter autocommands
134.914  000.055  000.055: sourcing /usr/share/nvim/runtime/ftplugin/conf.vim
136.306  000.068  000.068: sourcing /usr/share/nvim/runtime/syntax/conf.vim
145.236  014.569: editing files in windows
145.885  000.415  000.415: require('alpha')
145.943  000.055  000.055: require('alpha.themes.dashboard')
156.866  011.160: VimEnter autocommands
156.913  000.047: UIEnter autocommands
156.916  000.003: before starting main loop
159.794  000.258  000.258: require('bufferline.tabpages')
160.073  000.238  000.238: require('bufferline.models')
160.160  000.085  000.085: require('bufferline.pick')
160.299  000.138  000.138: require('bufferline.duplicates')
160.538  000.228  000.228: require('bufferline.diagnostics')
160.929  000.114  000.114: require('bufferline.numbers')
161.116  000.109  000.109: require('bufferline.sorters')
161.257  000.103  000.103: require('bufferline.offset')
161.348  000.088  000.088: require('bufferline.custom_area')
209.730  051.454: first screen update
209.735  000.005: --- NVIM STARTED ---

--- Startup times for process: Primary/TUI ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.167  000.165: event init
000.237  000.070: early init
000.272  000.036: locale set
000.308  000.036: init first window
000.524  000.216: inits 1
000.530  000.006: window checked
000.589  000.059: parsing arguments
000.951  000.029  000.029: require('vim.shared')
001.012  000.030  000.030: require('vim.inspect')
001.063  000.041  000.041: require('vim._options')
001.064  000.111  000.040: require('vim._editor')
001.065  000.160  000.020: require('vim._init_packages')
001.066  000.317: init lua interpreter
001.541  000.475: --- NVIM STARTED ---

--- Startup times for process: Embedded ---

times in msec
 clock   self+sourced   self:  sourced script
 clock   elapsed:              other lines

000.002  000.002: --- NVIM STARTING ---
000.128  000.126: event init
000.177  000.049: early init
000.206  000.028: locale set
000.232  000.027: init first window
000.513  000.281: inits 1
000.529  000.016: window checked
000.619  000.090: parsing arguments
001.235  000.050  000.050: require('vim.shared')
001.331  000.039  000.039: require('vim.inspect')
001.366  000.025  000.025: require('vim._options')
001.367  000.128  000.064: require('vim._editor')
001.367  000.233  000.054: require('vim._init_packages')
001.368  000.516: init lua interpreter
001.436  000.068: expanding arguments
001.454  000.017: inits 2
001.625  000.172: init highlight
001.626  000.001: waiting for UI
001.743  000.117: done waiting for UI
001.746  000.003: clear screen
001.793  000.003  000.003: require('vim.keymap')
002.163  000.416  000.413: require('vim._defaults')
002.165  000.003: init default mappings & autocommands
002.397  000.034  000.034: sourcing /usr/share/nvim/runtime/ftplugin.vim
002.439  000.015  000.015: sourcing /usr/share/nvim/runtime/indent.vim
002.516  000.047  000.047: sourcing /usr/share/nvim/archlinux.lua
002.519  000.065  000.018: sourcing /etc/xdg/nvim/sysinit.vim
003.150  000.415  000.415: require('KaGaSi.core.keymap')
004.215  001.063  001.063: require('KaGaSi.core.options')
004.288  000.071  000.071: require('KaGaSi.core.fold')
004.324  000.034  000.034: require('KaGaSi.core.filepos')
004.348  000.023  000.023: require('KaGaSi.core.pdf')
004.349  001.802  000.196: require('KaGaSi.core')
004.679  000.259  000.259: require('lazy')
004.712  000.017  000.017: require('ffi')
004.763  000.029  000.029: require('vim.fs')
004.908  000.142  000.142: require('vim.uri')
004.915  000.202  000.030: require('vim.loader')
005.079  000.140  000.140: require('lazy.stats')
005.168  000.069  000.069: require('lazy.core.util')
005.243  000.073  000.073: require('lazy.core.config')
005.390  000.064  000.064: require('lazy.core.handler')
005.497  000.106  000.106: require('lazy.core.plugin')
005.510  000.265  000.095: require('lazy.core.loader')
006.074  000.095  000.095: require('KaGaSi.plugins.alpha')
006.177  000.060  000.060: require('KaGaSi.plugins.auto-session')
006.272  000.088  000.088: require('KaGaSi.plugins.autopairs')
006.346  000.064  000.064: require('KaGaSi.plugins.bigfile')
006.412  000.057  000.057: require('KaGaSi.plugins.bufferline')
006.469  000.049  000.049: require('KaGaSi.plugins.colorscheme')
006.504  000.029  000.029: require('KaGaSi.plugins.comment')
006.538  000.024  000.024: require('KaGaSi.plugins.dressing')
006.577  000.033  000.033: require('KaGaSi.plugins.gitsigns')
006.608  000.026  000.026: require('KaGaSi.plugins.indent-blankline')
006.681  000.069  000.069: require('KaGaSi.plugins.leap')
006.758  000.058  000.058: require('KaGaSi.plugins.leap-spooky')
006.850  000.065  000.065: require('KaGaSi.plugins.lspconfig')
006.891  000.025  000.025: require('KaGaSi.plugins.ltex_extra')
006.919  000.026  000.026: require('KaGaSi.plugins.lualine')
006.956  000.029  000.029: require('KaGaSi.plugins.mason')
006.989  000.025  000.025: require('KaGaSi.plugins.numbertoggle')
007.036  000.043  000.043: require('KaGaSi.plugins.nvim-cmp')
007.380  000.046  000.046: require('KaGaSi.plugins.nvim-tree')
007.469  000.066  000.066: require('KaGaSi.plugins.surround')
007.528  000.048  000.048: require('KaGaSi.plugins.telescope')
007.571  000.030  000.030: require('KaGaSi.plugins.todo-comments')
007.605  000.027  000.027: require('KaGaSi.plugins.toggleterm')
007.643  000.034  000.034: require('KaGaSi.plugins.treesitter')
007.688  000.032  000.032: require('KaGaSi.plugins.trouble')
007.726  000.026  000.026: require('KaGaSi.plugins.vim-maximizer')
007.776  000.033  000.033: require('KaGaSi.plugins.vimtex')
008.797  000.974  000.974: require('vim.filetype')
008.800  001.017  000.043: require('KaGaSi.plugins.vimwiki')
008.846  000.030  000.030: require('KaGaSi.plugins.which_key')
009.042  000.047  000.047: require('lazy.core.handler.event')
009.067  000.023  000.023: require('lazy.core.handler.ft')
009.093  000.025  000.025: require('lazy.core.handler.keys')
009.116  000.021  000.021: require('lazy.core.handler.cmd')
009.497  000.016  000.016: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/cls.vim
009.526  000.010  000.010: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tex.vim
009.549  000.008  000.008: sourcing /home/gary/.local/share/nvim/lazy/vimtex/ftdetect/tikz.vim
009.951  000.122  000.122: sourcing /usr/share/nvim/runtime/filetype.lua
010.180  000.105  000.105: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/scroll.vim
010.244  000.044  000.044: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/plugin/themer.lua
011.303  000.363  000.363: require('themer')
011.346  000.040  000.040: require('themer.config')
011.397  000.033  000.033: require('themer.modules.core')
011.422  000.024  000.024: require('themer.modules.core.api')
011.470  000.047  000.047: require('themer.modules.themes.sonokai_deep')
011.507  000.036  000.036: require('themer.modules.core.utils')
011.606  000.097  000.097: require('themer.modules.core.mapper')
012.977  002.065  001.425: sourcing /home/gary/.local/share/nvim/lazy/themer.lua/colors/themer_sonokai_deep.lua
013.374  000.226  000.226: sourcing /home/gary/.local/share/nvim/lazy/leap.nvim/plugin/init.lua
013.906  000.457  000.457: require('leap-spooky')
015.394  000.026  000.026: sourcing /home/gary/.local/share/nvim/lazy/nvim-web-devicons/plugin/nvim-web-devicons.vim
015.727  000.025  000.025: require('bufferline.lazy')
015.765  000.034  000.034: require('bufferline.commands')
015.824  000.057  000.057: require('bufferline.config')
015.827  000.362  000.246: require('bufferline')
015.931  000.040  000.040: require('bufferline.utils')
015.933  000.102  000.062: require('bufferline.groups')
015.966  000.026  000.026: require('bufferline.constants')
016.015  000.047  000.047: require('bufferline.colors')
016.218  000.039  000.039: require('bufferline.highlights')
016.261  000.003  000.003: require('vim.F')
016.629  000.029  000.029: require('vim.version')
017.339  000.774  000.744: require('bufferline.hover')
017.480  000.064  000.064: require('bufferline.ui')
017.774  000.054  000.054: require('lualine_require')
017.970  000.445  000.391: require('lualine')
018.001  000.029  000.029: require('lazy.status')
018.966  000.053  000.053: require('lualine.utils.mode')
023.667  000.053  000.053: require('nvim-tree.notify')
023.673  000.104  000.051: require('nvim-tree.events')
023.841  000.050  000.050: require('nvim-tree.log')
024.081  000.151  000.151: require('nvim-tree.iterators.node-iterator')
024.108  000.265  000.114: require('nvim-tree.utils')
024.197  000.088  000.088: require('nvim-tree.git.utils')
024.364  000.166  000.166: require('nvim-tree.git.runner')
024.610  000.245  000.245: require('nvim-tree.watcher')
024.682  000.070  000.070: require('nvim-tree.explorer.node')
024.689  000.956  000.072: require('nvim-tree.git')
024.750  000.060  000.060: require('nvim-tree.explorer.watch')
024.808  000.033  000.033: require('nvim-tree.explorer.node-builders')
024.835  000.026  000.026: require('nvim-tree.explorer.sorters')
024.886  000.024  000.024: require('nvim-tree.marks')
024.890  000.053  000.029: require('nvim-tree.explorer.filters')
024.982  000.065  000.065: require('nvim-tree.view')
024.985  000.095  000.030: require('nvim-tree.live-filter')
024.987  000.236  000.028: require('nvim-tree.explorer.explore')
025.018  000.030  000.030: require('nvim-tree.explorer.reload')
025.020  001.346  000.064: require('nvim-tree.explorer')
025.022  001.494  000.044: require('nvim-tree.core')
025.049  000.026  000.026: require('nvim-tree.renderer.components.padding')
025.080  000.031  000.031: require('nvim-tree.renderer.components.icons')
025.106  000.025  000.025: require('nvim-tree.renderer.components.full-name')
025.202  000.021  000.021: require('nvim-tree.enum')
025.248  000.044  000.044: require('nvim-tree.renderer.decorator')
025.250  000.104  000.039: require('nvim-tree.renderer.decorator.bookmarks')
025.300  000.025  000.025: require('nvim-tree.renderer.decorator.copied')
025.322  000.021  000.021: require('nvim-tree.renderer.decorator.cut')
025.369  000.024  000.024: require('nvim-tree.diagnostics')
025.468  000.097  000.097: require('vim.diagnostic')
025.471  000.148  000.027: require('nvim-tree.renderer.decorator.diagnostics')
025.506  000.034  000.034: require('nvim-tree.renderer.decorator.git')
025.553  000.024  000.024: require('nvim-tree.buffers')
025.555  000.048  000.023: require('nvim-tree.renderer.decorator.modified')
025.577  000.021  000.021: require('nvim-tree.renderer.decorator.opened')
025.580  000.473  000.071: require('nvim-tree.renderer.builder')
025.585  002.126  000.077: require('nvim-tree.renderer')
025.588  002.234  000.108: require('nvim-tree.lib')
025.642  000.054  000.054: require('nvim-tree.appearance')
025.762  000.021  000.021: require('nvim-tree.actions.finders.find-file')
025.786  000.023  000.023: require('nvim-tree.actions.finders.search-node')
025.787  000.067  000.023: require('nvim-tree.actions.finders')
025.863  000.022  000.022: require('nvim-tree.actions.reloaders')
025.866  000.057  000.035: require('nvim-tree.actions.fs.copy-paste')
025.890  000.023  000.023: require('nvim-tree.actions.fs.create-file')
025.913  000.023  000.023: require('nvim-tree.actions.fs.remove-file')
025.938  000.024  000.024: require('nvim-tree.actions.fs.rename-file')
025.963  000.024  000.024: require('nvim-tree.actions.fs.trash')
025.964  000.176  000.026: require('nvim-tree.actions.fs')
026.007  000.021  000.021: require('nvim-tree.actions.moves.item')
026.027  000.019  000.019: require('nvim-tree.actions.moves.parent')
026.047  000.019  000.019: require('nvim-tree.actions.moves.sibling')
026.048  000.083  000.024: require('nvim-tree.actions.moves')
026.091  000.023  000.023: require('nvim-tree.actions.node.file-popup')
026.125  000.033  000.033: require('nvim-tree.actions.node.open-file')
026.176  000.050  000.050: require('nvim-tree.actions.node.run-command')
026.199  000.022  000.022: require('nvim-tree.actions.node.system-open')
026.200  000.152  000.024: require('nvim-tree.actions.node')
026.246  000.024  000.024: require('nvim-tree.actions.root.change-dir')
026.266  000.019  000.019: require('nvim-tree.actions.root.dir-up')
026.267  000.066  000.023: require('nvim-tree.actions.root')
026.320  000.029  000.029: require('nvim-tree.actions.tree.find-file')
026.368  000.026  000.026: require('nvim-tree.actions.tree.modifiers.collapse-all')
026.397  000.028  000.028: require('nvim-tree.actions.tree.modifiers.expand-all')
026.424  000.026  000.026: require('nvim-tree.actions.tree.modifiers.toggles')
026.426  000.105  000.025: require('nvim-tree.actions.tree.modifiers')
026.451  000.025  000.025: require('nvim-tree.actions.tree.open')
026.477  000.025  000.025: require('nvim-tree.actions.tree.toggle')
026.479  000.211  000.028: require('nvim-tree.actions.tree')
026.480  000.782  000.028: require('nvim-tree.actions')
026.515  000.034  000.034: require('nvim-tree.appearance.diagnostics')
026.586  000.036  000.036: require('nvim-tree.keymap')
026.589  000.073  000.038: require('nvim-tree.help')
026.618  000.028  000.028: require('nvim-tree.marks.navigation')
026.645  000.026  000.026: require('nvim-tree.marks.bulk-delete')
026.689  000.043  000.043: require('nvim-tree.marks.bulk-trash')
026.723  000.034  000.034: require('nvim-tree.marks.bulk-move')
026.750  001.080  000.060: require('nvim-tree.api')
026.755  001.111  000.031: require('nvim-tree.commands')
026.785  000.029  000.029: require('nvim-tree.legacy')
026.807  004.094  000.665: require('nvim-tree')
028.098  000.491  000.491: require('nvim-web-devicons.icons-default')
028.281  000.767  000.275: require('nvim-web-devicons')
030.722  000.045  000.045: require('auto-session.logger')
030.731  000.147  000.102: require('auto-session.lib')
030.770  000.037  000.037: require('auto-session.autocmds')
030.794  000.546  000.362: require('auto-session')
031.188  000.024  000.024: sourcing /home/gary/.local/share/nvim/lazy/plenary.nvim/plugin/plenary.vim
031.446  000.011  000.011: sourcing /home/gary/.local/share/nvim/lazy/todo-comments.nvim/plugin/todo.vim
031.892  000.136  000.136: require('todo-comments.util')
031.917  000.217  000.081: require('todo-comments.config')
032.183  000.185  000.185: require('todo-comments.highlight')
032.200  000.282  000.097: require('todo-comments.jump')
032.201  000.744  000.246: require('todo-comments')
032.665  000.361  000.361: sourcing /home/gary/.local/share/nvim/lazy/telescope.nvim/plugin/telescope.lua
032.975  000.073  000.073: require('telescope._extensions')
032.979  000.286  000.213: require('telescope')
033.298  000.025  000.025: require('plenary.bit')
033.322  000.023  000.023: require('plenary.functional')
033.345  000.129  000.081: require('plenary.path')
033.356  000.166  000.037: require('plenary.strings')
033.415  000.058  000.058: require('telescope.deprecated')
033.550  000.066  000.066: require('plenary.log')
033.579  000.122  000.056: require('telescope.log')
033.688  000.023  000.023: require('plenary.compat')
033.693  000.062  000.040: require('plenary.job')
033.724  000.030  000.030: require('telescope.state')
033.733  000.152  000.060: require('telescope.utils')
033.736  000.320  000.045: require('telescope.sorters')
034.484  001.342  000.799: require('telescope.config')
034.561  000.031  000.031: require('plenary.window.border')
034.586  000.023  000.023: require('plenary.window')
034.607  000.021  000.021: require('plenary.popup.utils')
034.609  000.123  000.049: require('plenary.popup')
034.641  000.031  000.031: require('telescope.pickers.scroller')
034.671  000.029  000.029: require('telescope.actions.state')
034.702  000.030  000.030: require('telescope.actions.utils')
034.766  000.031  000.031: require('telescope.actions.mt')
034.773  000.070  000.039: require('telescope.actions.set')
034.834  000.031  000.031: require('telescope.config.resolve')
034.836  000.062  000.031: require('telescope.pickers.entry_display')
034.870  000.033  000.033: require('telescope.from_entry')
034.980  002.000  000.279: require('telescope.actions')
035.582  000.092  000.092: require('fzf_lib')
035.586  000.148  000.056: require('telescope._extensions.fzf')
035.659  004.846  001.272: require('telescope')
035.728  000.024  000.024: require('auto-session.session-lens.library')
035.753  000.024  000.024: require('auto-session.session-lens.actions')
035.779  000.118  000.070: require('auto-session.session-lens')
036.111  000.022  000.022: require('toggleterm.lazy')
036.133  000.019  000.019: require('toggleterm.constants')
036.187  000.053  000.053: require('toggleterm.terminal')
036.190  000.252  000.157: require('toggleterm')
036.244  000.021  000.021: require('toggleterm.colors')
036.273  000.028  000.028: require('toggleterm.utils')
036.276  000.085  000.036: require('toggleterm.config')
037.032  000.074  000.074: require('toggleterm.commandline')
037.131  000.023  000.023: sourcing /home/gary/.local/share/nvim/lazy/vim-numbertoggle/plugin/number_toggle.vim
038.305  000.057  000.057: require('mason-core.functional')
038.340  000.032  000.032: require('mason-core.path')
038.374  000.034  000.034: require('mason.settings')
038.388  000.176  000.054: require('mason-core.log')
038.389  000.209  000.032: require('mason-core.EventEmitter')
038.422  000.032  000.032: require('mason-core.optional')
038.530  000.047  000.047: require('mason-core.async')
038.557  000.025  000.025: require('mason-core.async.uv')
038.560  000.136  000.064: require('mason-core.fs')
038.597  000.036  000.036: require('mason-registry.sources')
038.658  000.026  000.026: require('mason-core.functional.data')
038.691  000.033  000.033: require('mason-core.functional.function')
038.703  000.104  000.045: require('mason-core.functional.list')
038.713  000.752  000.234: require('mason-registry')
038.716  001.041  000.289: require('mason-tool-installer')
038.726  001.082  000.041: sourcing /home/gary/.local/share/nvim/lazy/mason-tool-installer.nvim/plugin/mason-tool-installer.lua
038.880  000.028  000.028: require('mason-core.functional.relation')
038.942  000.032  000.032: require('mason-core.functional.logic')
038.950  000.151  000.091: require('mason-core.platform')
038.951  000.185  000.034: require('mason')
039.012  000.024  000.024: require('mason-lspconfig.settings')
039.015  000.063  000.039: require('mason-lspconfig')
039.265  000.048  000.048: require('mason-core.functional.string')
039.283  000.139  000.092: require('mason.api.command')
039.382  000.035  000.035: require('mason-lspconfig.notify')
039.384  000.083  000.048: require('mason-lspconfig.lspconfig_hook')
040.287  000.059  000.059: require('lsp-file-operations.log')
040.290  000.273  000.214: require('lsp-file-operations')
041.562  000.294  000.294: require('neodev')
041.611  000.046  000.046: require('neodev.config')
041.708  000.049  000.049: require('neodev.util')
041.710  000.092  000.043: require('neodev.lsp')
042.100  000.142  000.142: require('vim.lsp.log')
042.465  000.362  000.362: require('vim.lsp.protocol')
042.972  000.315  000.315: require('vim.lsp._snippet_grammar')
043.023  000.048  000.048: require('vim.highlight')
043.043  000.576  000.213: require('vim.lsp.util')
043.161  000.064  000.064: require('vim.lsp.sync')
043.168  000.124  000.060: require('vim.lsp._changetracking')
043.295  000.125  000.125: require('vim.lsp.rpc')
043.330  001.548  000.220: require('vim.lsp')
043.389  001.678  000.130: require('lspconfig.util')
043.821  000.250  000.250: sourcing /home/gary/.local/share/nvim/lazy/nvim-lspconfig/plugin/lspconfig.lua
044.217  000.164  000.164: require('lspconfig.async')
044.219  000.286  000.122: require('lspconfig.configs')
044.227  000.385  000.099: require('lspconfig')
044.567  000.236  000.236: require('cmp_nvim_lsp.source')
044.569  000.341  000.106: require('cmp_nvim_lsp')
045.008  000.198  000.198: require('mason-core.functional.table')
045.035  000.410  000.212: require('mason-lspconfig.mappings.server')
045.219  000.093  000.093: require('lspconfig.server_configurations.cssls')
045.507  000.043  000.043: require('lspconfig.manager')
045.563  000.051  000.051: require('lspconfig.server_configurations.html')
045.736  000.034  000.034: require('lspconfig.server_configurations.clangd')
045.936  000.032  000.032: require('lspconfig.server_configurations.bashls')
046.075  000.034  000.034: require('lspconfig.server_configurations.pyright')
046.214  000.039  000.039: require('lspconfig.server_configurations.vimls')
046.357  000.036  000.036: require('lspconfig.server_configurations.lua_ls')
046.512  000.081  000.081: require('lspconfig.server_configurations.ltex')
046.602  007.217  003.005: require('lspconfig.util')
046.633  000.029  000.029: require('mason-lspconfig.server_config_extensions')
046.670  000.035  000.035: require('lspconfig.server_configurations.omnisharp')
046.724  000.028  000.028: require('mason-lspconfig.ensure_installed')
046.798  000.028  000.028: require('mason-core.result')
046.896  000.046  000.046: require('mason-core.process')
046.937  000.040  000.040: require('mason-core.spawn')
046.939  000.114  000.029: require('mason-core.managers.powershell')
046.961  000.021  000.021: require('mason.version')
046.962  000.162  000.027: require('mason-core.fetch')
046.987  000.024  000.024: require('mason-core.providers')
047.137  000.075  000.075: require('mason-core.purl')
047.152  000.140  000.065: require('mason-core.package')
047.295  000.045  000.045: require('mason-core.installer.registry.expr')
047.301  000.097  000.052: require('mason-core.installer.registry.link')
047.540  000.028  000.028: require('mason-core.receipt')
047.552  000.085  000.057: require('mason-core.installer.context')
047.581  000.028  000.028: require('mason-core.async.control')
047.609  000.027  000.027: require('mason-core.installer.linker')
047.611  000.193  000.053: require('mason-core.installer')
047.620  000.250  000.057: require('mason-core.installer.managers.std')
047.621  000.319  000.069: require('mason-core.installer.registry.schemas')
047.647  000.025  000.025: require('mason-core.installer.registry.util')
047.652  000.498  000.057: require('mason-core.installer.registry')
047.653  000.665  000.028: require('mason-registry.sources.util')
047.658  000.918  000.039: require('mason-registry.sources.github')
050.553  000.025  000.025: require('mason-core.functional.number')
050.569  000.082  000.057: require('mason-lspconfig.api.command')
051.111  000.037  000.037: sourcing /usr/share/nvim/runtime/plugin/editorconfig.lua
051.231  000.100  000.100: sourcing /usr/share/nvim/runtime/plugin/gzip.vim
051.283  000.037  000.037: sourcing /usr/share/nvim/runtime/plugin/man.lua
051.708  000.136  000.136: sourcing /usr/share/nvim/runtime/pack/dist/opt/matchit/plugin/matchit.vim
051.718  000.419  000.283: sourcing /usr/share/nvim/runtime/plugin/matchit.vim
051.964  000.229  000.229: sourcing /usr/share/nvim/runtime/plugin/matchparen.vim
052.000  000.009  000.009: sourcing /usr/share/nvim/runtime/plugin/netrwPlugin.vim
052.103  000.088  000.088: sourcing /usr/share/nvim/runtime/plugin/osc52.lua
052.210  000.090  000.090: sourcing /usr/share/nvim/runtime/plugin/rplugin.vim
052.272  000.044  000.044: sourcing /usr/share/nvim/runtime/plugin/shada.vim
052.298  000.009  000.009: sourcing /usr/share/nvim/runtime/plugin/spellfile.vim
052.363  000.050  000.050: sourcing /usr/share/nvim/runtime/plugin/tarPlugin.vim
052.424  000.042  000.042: sourcing /usr/share/nvim/runtime/plugin/tohtml.lua
052.453  000.011  000.011: sourcing /usr/share/nvim/runtime/plugin/tutor.vim
052.541  000.068  000.068: sourcing /usr/share/nvim/runtime/plugin/zipPlugin.vim
052.707  000.031  000.031: sourcing /home/gary/.local/share/nvim/lazy/cmp-nvim-lsp/after/plugin/cmp_nvim_lsp.lua
052.918  000.090  000.090: require('bigfile.features')
052.922  000.152  000.062: require('bigfile')
052.936  000.194  000.042: sourcing /home/gary/.local/share/nvim/lazy/bigfile.nvim/after/plugin/bigfile.lua
052.970  048.620  018.081: require('KaGaSi.lazy')
052.971  050.433  000.011: sourcing /home/gary/.config/nvim/init.lua
052.975  000.264: sourcing vimrc file(s)
053.302  000.044  000.044: sourcing /usr/share/nvim/runtime/filetype.lua
053.479  000.043  000.043: sourcing /usr/share/nvim/runtime/syntax/synload.vim
054.256  000.905  000.862: sourcing /usr/share/nvim/runtime/syntax/syntax.vim
054.276  000.352: inits 3
056.046  001.771: reading ShaDa
056.070  000.023: making windows
056.682  000.026  000.026: require('nvim-surround.input')
056.706  000.022  000.022: require('nvim-surround.functional')
056.712  000.137  000.089: require('nvim-surround.config')
056.716  000.195  000.058: require('nvim-surround.buffer')
056.738  000.022  000.022: require('nvim-surround.cache')
056.767  000.027  000.027: require('nvim-surround.utils')
056.770  000.473  000.229: require('nvim-surround')
057.440  000.035  000.035: require('gitsigns.async')
057.472  000.029  000.029: require('gitsigns.debug.log')
057.556  000.083  000.083: require('gitsigns.config')
057.559  000.477  000.328: require('gitsigns')
057.638  000.053  000.053: require('gitsigns.highlight')
058.086  000.088  000.088: require('gitsigns.util')
058.114  000.026  000.026: require('gitsigns.system')
058.140  000.024  000.024: require('gitsigns.message')
058.168  000.028  000.028: require('gitsigns.git.version')
058.179  000.260  000.094: require('gitsigns.git')
058.252  000.026  000.026: require('gitsigns.cache')
058.279  000.025  000.025: require('gitsigns.signs')
058.304  000.024  000.024: require('gitsigns.status')
058.326  000.021  000.021: require('gitsigns.debounce')
058.348  000.021  000.021: require('gitsigns.diff')
058.378  000.029  000.029: require('gitsigns.hunks')
058.383  000.204  000.057: require('gitsigns.manager')
058.387  000.520  000.056: require('gitsigns.attach')
058.439  000.038  000.038: require('gitsigns.current_line_blame')
058.529  000.043  000.043: require('vim._system')
061.927  000.041  000.041: require('nvim-treesitter.utils')
062.355  000.029  000.029: require('vim.treesitter.language')
062.391  000.027  000.027: require('vim.func')
062.418  000.025  000.025: require('vim.func._memoize')
062.440  000.171  000.090: require('vim.treesitter.query')
062.490  000.033  000.033: require('vim.treesitter._range')
062.499  000.297  000.093: require('vim.treesitter.languagetree')
062.503  000.345  000.048: require('vim.treesitter')
063.085  001.155  000.810: require('nvim-treesitter.parsers')
063.511  000.065  000.065: require('nvim-treesitter.compat')
063.671  000.111  000.111: require('nvim-treesitter.ts_utils')
063.680  000.168  000.057: require('nvim-treesitter.tsrange')
063.728  000.047  000.047: require('nvim-treesitter.caching')
063.740  000.478  000.199: require('nvim-treesitter.query')
063.759  000.598  000.120: require('nvim-treesitter.configs')
063.766  000.679  000.081: require('nvim-treesitter.info')
063.865  000.099  000.099: require('nvim-treesitter.shell_command_selectors')
063.912  003.577  001.602: require('nvim-treesitter.install')
063.956  000.042  000.042: require('nvim-treesitter.statusline')
064.031  000.074  000.074: require('nvim-treesitter.query_predicates')
064.032  003.943  000.250: require('nvim-treesitter')
064.303  004.258  000.315: sourcing /home/gary/.local/share/nvim/lazy/nvim-treesitter/plugin/nvim-treesitter.lua
065.835  000.077  000.077: require('nvim-treesitter.indent')
066.225  000.082  000.082: require('nvim-treesitter.locals')
066.230  000.139  000.057: require('nvim-treesitter.incremental_selection')
066.383  000.055  000.055: require('nvim-treesitter.highlight')
067.165  000.127  000.127: require('ibl.utils')
067.171  000.210  000.083: require('ibl.config')
067.285  000.061  000.061: require('ibl.indent')
067.298  000.126  000.065: require('ibl.hooks')
067.301  000.421  000.086: require('ibl.highlights')
067.346  000.045  000.045: require('ibl.autocmds')
067.399  000.051  000.051: require('ibl.inlay_hints')
067.463  000.063  000.063: require('ibl.virt_text')
067.782  000.240  000.240: require('ibl.scope_languages')
067.786  000.322  000.081: require('ibl.scope')
067.793  001.232  000.329: require('ibl')
067.813  001.301  000.069: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/plugin/commands.lua
068.009  000.126  000.126: require('vim.iter')
068.265  000.214  000.214: require('vim.lsp.handlers')
069.501  000.100  000.100: sourcing /home/gary/.local/share/nvim/lazy/nvim-ts-context-commentstring/plugin/ts_context_commentstring.lua
070.118  000.049  000.049: require('Comment.config')
070.384  000.202  000.202: require('Comment.ft')
070.395  000.275  000.072: require('Comment.utils')
070.523  000.127  000.127: require('Comment.opfunc')
070.646  000.121  000.121: require('Comment.extra')
070.650  001.020  000.448: require('Comment.api')
070.787  001.204  000.184: sourcing /home/gary/.local/share/nvim/lazy/Comment.nvim/plugin/Comment.lua
070.915  000.105  000.105: require('Comment')
071.029  000.107  000.107: require('ts_context_commentstring.integrations.comment_nvim')
073.766  000.166  000.166: sourcing /usr/share/nvim/runtime/autoload/provider/clipboard.vim
076.924  001.988  001.988: require('vim.filetype.detect')
077.534  000.053  000.053: sourcing /usr/share/nvim/runtime/ftplugin/c.vim
077.674  000.126  000.126: sourcing /usr/share/nvim/runtime/ftplugin/c.lua
078.075  000.168  000.168: sourcing /home/gary/.local/share/nvim/lazy/indent-blankline.nvim/after/ftplugin/c.lua
078.363  000.012  000.012: sourcing /usr/share/nvim/runtime/indent/c.vim
079.361  000.419  000.419: require('vim.lsp.client')
080.003  000.224  000.224: require('vim.glob')
080.039  000.450  000.226: require('vim.lsp._dynamic')
089.404  000.620  000.620: sourcing /usr/share/nvim/runtime/syntax/c.vim
117.286  000.154  000.154: require('vim.treesitter.highlighter')
124.171  000.033  000.033: require('ts_context_commentstring.utils')
124.220  000.045  000.045: require('ts_context_commentstring.config')
124.222  000.154  000.076: require('ts_context_commentstring.internal')
124.386  000.062  000.062: require('editorconfig')
124.478  054.747: opening buffers
124.525  000.037  000.037: require('bufferline.state')
126.376  001.861: BufEnter autocommands
130.318  000.029  000.029: sourcing /usr/share/nvim/runtime/ftplugin/conf.vim
131.517  000.046  000.046: sourcing /usr/share/nvim/runtime/syntax/conf.vim
139.436  012.986: editing files in windows
140.092  000.484  000.484: require('alpha')
140.161  000.067  000.067: require('alpha.themes.dashboard')
153.739  013.752: VimEnter autocommands
153.849  000.110: UIEnter autocommands
153.859  000.010: before starting main loop
157.834  000.512  000.512: require('bufferline.tabpages')
158.237  000.370  000.370: require('bufferline.models')
158.404  000.165  000.165: require('bufferline.pick')
158.564  000.158  000.158: require('bufferline.duplicates')
158.867  000.290  000.290: require('bufferline.diagnostics')
159.573  000.187  000.187: require('bufferline.numbers')
159.951  000.206  000.206: require('bufferline.sorters')
160.217  000.223  000.223: require('bufferline.offset')
160.385  000.165  000.165: require('bufferline.custom_area')
201.880  045.744: first screen update
201.886  000.006: --- NVIM STARTED ---

