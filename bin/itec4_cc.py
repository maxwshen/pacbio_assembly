# Combines contigs
import itec4

def main():
  contigs_fold = sys.argv[1]
  itec4.combine_contigs(contigs_fold)

if __name__ == '__main__':
  start = datetime.datetime.now()
  print 'Start:', start, '\n'
  main()
  end = datetime.datetime.now()
  print '\n\nEnd:', end, '\nTotal:', end - start