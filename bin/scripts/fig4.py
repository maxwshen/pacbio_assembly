inp_file = '/home/mshen/research/figure_data/raw_3deg_fig4.txt'
out_file = '/home/mshen/research/figure_data/3deg_fig4_4.gv'

creads = []
with open(inp_file) as f:
  for i, line in enumerate(f):
    items = line.strip('[]')
    items = items.replace('\'', '').split(', ')
    creads.append(items[1:-1])

with open(out_file, 'w') as f:
  f.write('digraph G {\n  rankdir=LR;\n')
  base_nodes = []
  for j in range(len(creads) / 8):
    cread = creads[j]
    for i in range(0, len(cread) - 1, 2):
      if j == 0:
        base_nodes += [cread[i], cread[i+2]]
        edge = '  ' + cread[i] + ' -> ' + cread[i+2] + '[color=red,label=\"' + cread[i+1] + '\"];\n'
      else:
        edge = '  ' + cread[i] + ' -> ' + cread[i+2] + '[label=\"' + cread[i+1] + '\"];\n'
      f.write(edge)
  # f.write('  { rank=same, ' + ', '.join(base_nodes) + ' }\n')
  f.write('}')
