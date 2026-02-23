import groovy.transform.CompileDynamic

import java.nio.file.Path
import java.text.SimpleDateFormat

import groovy.yaml.YamlSlurper

/**
  * A series of utility functions for testing the pipeline.
  */
@CompileDynamic
class Utils {

    static Object getRecursiveFileNames(Path fileOrDir, CharSequence outputDir) {
        if (file(fileOrDir.toString()).directory) {
            return fileOrDir.list().collect { Path f -> getRecursiveFileNames(f, outputDir) }
        }
        return fileOrDir.toString().replace("${outputDir}/", '')
    }

    static String getDynamicOutputName() {
        Map nfcoreYaml = new YamlSlurper().parseText(file('.nf-core.yml').text) as Map
        String date = new SimpleDateFormat('yyyy_MM_dd', Locale.FRANCE).format(new Date()) as String
        String version = ((nfcoreYaml.get('template') as Map).get('version') as String).replace('.', '_')
        return "v${version}_${date}" as String
    }

}
