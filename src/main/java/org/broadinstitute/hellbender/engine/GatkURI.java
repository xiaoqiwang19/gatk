package org.broadinstitute.hellbender.engine;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import org.broadinstitute.barclay.argparser.TaggedArgument;
import org.broadinstitute.barclay.argparser.TaggedArgumentParser;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;

import java.nio.file.Path;
import java.util.Map;

// TODO: need a way to tag GCS-enabled args
// TODO: should we do anything with query params (i.e. -> tags) ?
public class GatkURI extends HtsSimpleURI implements TaggedArgument {

    private String tagName;
    private Map<String, String> tagAttributes;

    public GatkURI(final String uriString) {
        super(uriString);
    }

    @Override
    public Path toPath() {
        // special case GCS, in case the filesystem provider wasn't installed properly but is available.
        if (CloudStorageFileSystem.URI_SCHEME.equals(getURI().getScheme())) {
            final Path tempPath = BucketUtils.getPathOnGcs(getURIString());
            setCachedPath(tempPath);
            return tempPath;
        } else {
            return super.toPath();
        }
    }

    @Override
    public void setTag(String tagName) {
        this.tagName = tagName;
    }

    @Override
    public String getTag() {
        return tagName;
    }

    @Override
    public void setTagAttributes(final Map<String, String> attributes) {
            this.tagAttributes = attributes;
    }

    @Override
    public Map<String, String> getTagAttributes() {
        return tagAttributes;
    }

    @Override
    public String toString() {
        if (getTagAttributes() != null) {
            return super.toString() + TaggedArgumentParser.getDisplayString("", this);
        } else {
            return super.toString();
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof GatkURI)) return false;
        if (!super.equals(o)) return false;

        GatkURI gatkURI = (GatkURI) o;

        if (!tagName.equals(gatkURI.tagName)) return false;
        return getTagAttributes().equals(gatkURI.getTagAttributes());
    }

    @Override
    public int hashCode() {
        int result = super.hashCode();
        result = 31 * result + tagName.hashCode();
        result = 31 * result + getTagAttributes().hashCode();
        return result;
    }
}
